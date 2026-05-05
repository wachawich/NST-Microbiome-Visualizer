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
        จัดการชื่อ Taxonomy ตามเงื่อนไข:
        1. รวมกลุ่มข้อมูลตาม Rank ที่เลือก
        2. ถ้าชื่อเป็น unknown/uncultured ให้ใช้ Rank สูงกว่าแทน
        """
        target_rank = target_rank.lower()
        rank_map = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
        full_rank_map = {'domain': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6}
        
        if target_rank in full_rank_map:
            rank_idx = full_rank_map[target_rank]
        else:
            rank_idx = rank_map.get(target_rank[0], 5)
            
        taxa_mapping = {}
        for col in self.taxa_cols:
            parts = col.split(';')
            name = 'Unknown'
            for i in range(min(rank_idx, len(parts) - 1), -1, -1):
                part = parts[i].strip()
                val = part[3:] if len(part) > 3 and part[1:3] == '__' else part
                if val.lower() not in ['', 'incertae_sedis', 'unknown', 'uncultured', 'unassigned'] and not part.endswith('__'):
                    name = part
                    break
            taxa_mapping[col] = name
            
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

    def heatmap(self, taxonomic_rank='genus', scaling='z-score', 
                       clustering=True, file_format='png', filename='heatmap'):
        """ฟังก์ชันสำหรับสร้าง Heatmap"""
        return plt
        
    def export_heatmap(self, taxonomic_rank='genus', scaling='z-score', 
                       clustering=True, file_format='png', filename='heatmap'):
        """ฟังก์ชันสำหรับเซฟ Heatmap"""
        plt.savefig(f"{filename}.{file_format}")
        plt.close()
        
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