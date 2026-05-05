import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class MicrobiomeVisualizer:
    def __init__(self, file_path):
        """โหลดข้อมูลและเตรียมโครงสร้างเบื้องต้น"""
        self.df = pd.read_csv(file_path)
        if 'index' in self.df.columns:
            self.df = self.df.set_index('index')
        self.metadata_col = self.df.columns[-1]  # คอลัมน์ Organic/Chemical
        self.taxa_cols = self.df.columns[:-1].tolist()
        
    def _parse_taxonomy(self, target_rank):
        """
        จัดกลุ่มตาม Rank ที่เลือก:
        - ถ้าระดับนั้นมีชื่อใช้งานได้ → ใช้ชื่อนั้น
        - ถ้าระบุไม่ได้ (Incertae_Sedis / uncultured / unknown / __ / ฯลฯ)
          → ยุบรวมเป็น 'Unclassified' (ไม่ fallback ไประดับสูงกว่า)
        """
        target_rank = target_rank.lower()
        rank_map = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
        full_rank_map = {'domain': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6}

        if target_rank in full_rank_map:
            rank_idx = full_rank_map[target_rank]
        else:
            rank_idx = rank_map.get(target_rank[0], 5)

        invalid = {'', 'incertae_sedis', 'unknown', 'uncultured', 'unassigned'}
        taxa_mapping = {}
        for col in self.taxa_cols:
            parts = col.split(';')
            if rank_idx >= len(parts):
                taxa_mapping[col] = 'Unclassified'
                continue
            part = parts[rank_idx].strip()
            val = part[3:] if len(part) > 3 and part[1:3] == '__' else part
            if val.lower() in invalid or part.endswith('__'):
                taxa_mapping[col] = 'Unclassified'
            else:
                taxa_mapping[col] = part

        return taxa_mapping

    def _apply_cutoff(self, df_relative, cutoff_percent):
        """กรองจุลินทรีย์ที่ต่ำกว่า cutoff และยุบรวมเป็น 'Others'"""
        keep_taxa = df_relative.columns[df_relative.max() >= cutoff_percent]
        
        df_filtered = df_relative[keep_taxa].copy()
        df_filtered['Others'] = df_relative.drop(columns=keep_taxa).sum(axis=1)
        
        sorted_taxa = df_filtered.drop(columns=['Others']).mean().sort_values(ascending=False).index.tolist()
        if 'Others' in df_filtered.columns and df_filtered['Others'].sum() > 0:
            sorted_taxa.append('Others')
        else:
            df_filtered = df_filtered.drop(columns=['Others'], errors='ignore')
            
        return df_filtered[sorted_taxa]
        
    # Function นี้ของ stacked bar plot
    def get_relative_abundance(self, taxonomic_rank='genus', percent_cutoff=1.0):
        """คำนวณและเตรียมข้อมูลสำหรับ plot"""

        # Map columns to new taxonomy names
        taxa_map = self._parse_taxonomy(taxonomic_rank)
        
        # Calculate relative abundance per sample
        df_counts = self.df[self.taxa_cols]
        df_rel = df_counts.div(df_counts.sum(axis=1), axis=0) * 100
        
        # Rename columns and group by new taxonomy
        df_rel = df_rel.rename(columns=taxa_map)
        df_rel = df_rel.T.groupby(df_rel.columns).sum().T
        
        # Group by metadata (Farm type)
        df_rel['Group'] = self.df[self.metadata_col]
        df_grouped = df_rel.groupby('Group').mean()
        
        # Apply Cutoff
        df_final = self._apply_cutoff(df_grouped, percent_cutoff)
        return df_final
    
    def stacked_bar(self, taxonomic_rank='genus', percent_cutoff=1.0, 
                           file_format='pdf', filename='stacked_bar'):
        """ฟังก์ชันสำหรับสร้าง Stacked Bar Plot"""

        df_final = self.get_relative_abundance(taxonomic_rank, percent_cutoff)
        ax = df_final.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20', edgecolor='white', linewidth=0.5)
        
        plt.title(f'Microbiome Composition ({taxonomic_rank.capitalize()} level) - {percent_cutoff}% Cutoff', pad=20)
        plt.ylabel('Relative Abundance (%)')
        plt.xlabel('Farm Type')
        plt.xticks(rotation=0)
        plt.legend(title='Taxa', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        return plt

    def export_stacked_bar(self, taxonomic_rank='genus', percent_cutoff=1.0, 
                           file_format='pdf', filename='stacked_bar'):
        """ฟังก์ชันสำหรับเซฟ Stacked Bar Plot"""
        plt = self.stacked_bar(taxonomic_rank, percent_cutoff)
        plt.savefig(f"{filename}.{file_format}", bbox_inches='tight')
        plt.close()

    def get_abundance_matrix(self, taxonomic_rank='genus', percent_cutoff=1.0):
        """เตรียมเมทริกซ์ Relative Abundance ระดับ sample (rows = samples, cols = taxa)"""
        taxa_map = self._parse_taxonomy(taxonomic_rank)

        df_counts = self.df[self.taxa_cols]
        df_rel = df_counts.div(df_counts.sum(axis=1), axis=0) * 100

        df_rel = df_rel.rename(columns=taxa_map)
        df_rel = df_rel.T.groupby(df_rel.columns).sum().T

        if percent_cutoff > 0:
            keep_taxa = df_rel.columns[df_rel.max() >= percent_cutoff]
            df_rel = df_rel[keep_taxa]

        sorted_taxa = df_rel.mean().sort_values(ascending=False).index.tolist()
        return df_rel[sorted_taxa]

    def heatmap(self, taxonomic_rank='genus', percent_cutoff=1.0, scaling='z-score',
                clustering=True):
        """ฟังก์ชันสำหรับสร้าง Heatmap
        - แกน X = Samples (Farm ID), แกน Y = Taxa
        - scaling: 'z-score', 'log', หรือ None (ใช้ค่า relative abundance ดิบ)
        - clustering=True ใช้ seaborn.clustermap เพื่อแสดง dendrogram และจัดกลุ่มแบบ hierarchical
        """
        df = self.get_abundance_matrix(taxonomic_rank, percent_cutoff)

        if scaling == 'log':
            data = np.log1p(df)
            cbar_label = 'log(1 + Relative Abundance %)'
            cmap = 'viridis'
            center = None
        elif scaling == 'z-score':
            std = df.std(axis=0).replace(0, 1)
            data = (df - df.mean(axis=0)) / std
            cbar_label = 'Z-score (per Taxon)'
            cmap = 'vlag'
            center = 0
        else:
            data = df.copy()
            cbar_label = 'Relative Abundance (%)'
            cmap = 'viridis'
            center = None

        # transpose ให้ rows = taxa, cols = samples ตามสเปค
        data_T = data.T

        groups = self.df[self.metadata_col].loc[data.index]
        unique_groups = sorted(groups.unique())
        default_palette = {'Organic': '#2ca02c', 'Chemical': '#d62728'}
        group_palette = {g: default_palette.get(g, c)
                         for g, c in zip(unique_groups,
                                         sns.color_palette('Set2', len(unique_groups)))}
        col_colors = groups.map(group_palette)

        title = (f'Microbiome Heatmap ({taxonomic_rank.capitalize()} level, '
                 f'{percent_cutoff}% cutoff, scaling={scaling})')

        plt.close('all')

        if clustering:
            g = sns.clustermap(
                data_T,
                cmap=cmap, center=center,
                col_colors=col_colors,
                figsize=(max(8, 0.6 * data_T.shape[1] + 4),
                         max(6, 0.25 * data_T.shape[0] + 2)),
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': cbar_label},
                linewidths=0,
            )
            g.ax_heatmap.set_xlabel('Sample')
            g.ax_heatmap.set_ylabel('Taxa')
            import matplotlib.patches as mpatches
            patches = [mpatches.Patch(color=group_palette[name], label=name)
                       for name in unique_groups]
            g.figure.subplots_adjust(top=0.90)
            g.figure.suptitle(title, y=0.98, fontsize=11)
            g.figure.legend(handles=patches, title=self.metadata_col,
                            loc='upper center', bbox_to_anchor=(0.5, 0.94),
                            ncol=len(group_palette), frameon=False)
        else:
            _, ax = plt.subplots(
                figsize=(max(8, 0.6 * data_T.shape[1] + 4),
                         max(6, 0.25 * data_T.shape[0] + 2)))
            sns.heatmap(data_T, cmap=cmap, center=center, ax=ax,
                        cbar_kws={'label': cbar_label},
                        xticklabels=True, yticklabels=True)
            ax.set_xlabel('Sample')
            ax.set_ylabel('Taxa')
            ax.set_title(title)
            plt.tight_layout()

        return plt

    def export_heatmap(self, taxonomic_rank='genus', percent_cutoff=1.0, scaling='z-score',
                       clustering=True, file_format='png', filename='heatmap'):
        """ฟังก์ชันสำหรับเซฟ Heatmap"""
        plt = self.heatmap(taxonomic_rank, percent_cutoff, scaling, clustering)
        plt.savefig(f"{filename}.{file_format}", bbox_inches='tight')
        plt.close('all')
        
    def bubble_plot(self, taxonomic_rank='genus', percent_cutoff=2.0,
                           file_format='tiff', filename='bubble'):
        """ฟังก์ชันสำหรับสร้าง Bubble Plot"""
        return plt

    def export_bubble_plot(self, taxonomic_rank='genus', percent_cutoff=2.0,
                           file_format='tiff', filename='bubble'):
        """ฟังก์ชันสำหรับเซฟ Bubble Plot"""
        plt = self.bubble_plot()
        plt.savefig(f"{filename}.{file_format}")
        plt.close()