# Microbiome Visualizer

Microbiome Visualizer เป็นคลาส Python ที่ออกแบบมาเพื่อช่วยในการวิเคราะห์ ประมวลผล และสร้างกราฟจากข้อมูลจุลินทรีย์ (Microbiome) โดยรองรับข้อมูลตาราง (Feature Table) ที่มีข้อมูล Taxonomy และข้อมูลกลุ่มตัวอย่าง (Metadata)

## ข้อมูล
- สูตรคำนวน Relative Abundance:

$$
\text{Relative Abundance} = \frac{\text{Abundance of a specific taxon}}{\text{Total abundance of all taxa in the sample}} \times 100
$$

- dataset : feature_table.csv


## วิธีการใช้งานเบื้องต้น (Quick Start)

### 1. โหลดข้อมูล
คุณสามารถเรียกใช้คลาสและนำเข้าข้อมูล CSV ได้ทันที (ข้อมูลในคอลัมน์แรกควรเป็น Index และคอลัมน์สุดท้ายควรเป็น Metadata ของตัวแปรต้น)
```python
from main import MicrobiomeVisualizer

# สร้างออบเจ็กต์ด้วยไฟล์ feature_table.csv
visualizer = MicrobiomeVisualizer('feature_table.csv')
```

### 2. การสร้างกราฟ Stacked Bar Plot
สามารถสร้าง Stacked Bar Plot เพื่อดูสัดส่วนสัมพัทธ์ (Relative Abundance) ของจุลินทรีย์ในระดับอนุกรมวิธาน (Taxonomic rank) ที่ต้องการได้
```python
import matplotlib.pyplot as plt

# ดูระดับ Genus และนำจุลินทรีย์ที่มีสัดส่วนเฉลี่ยน้อยกว่า 1% ไปรวมในกลุ่ม 'Others'
plt_fig = visualizer.stacked_bar(taxonomic_rank='genus', percent_cutoff=1.0)
plt_fig.show()
```

### 3. การบันทึกกราฟเป็นไฟล์ (Export)
รองรับการบันทึกกราฟเป็นไฟล์นามสกุลต่างๆ ได้ทันที เช่น `.pdf`, `.png`, `.tiff`, `.jpg`
```python
# บันทึกกราฟระดับ Phylum ตัดที่ 5% เป็นไฟล์ pdf
visualizer.export_stacked_bar(taxonomic_rank='phylum', percent_cutoff=5.0, 
                              file_format='pdf', filename='my_phylum_plot')
```

---

## การทำงานของฟังก์ชัน `_parse_taxonomy`

ฟังก์ชัน `_parse_taxonomy` เป็นฟังก์ชันหลักที่อยู่เบื้องหลังการจัดกลุ่มข้อมูลจุลินทรีย์ หน้าที่ของมันคือการ **แปลงสายอักขระ Taxonomy รูปแบบยาว ให้กลายเป็นชื่อที่เหมาะสมที่สุดในระดับชั้นที่เราระบุ** โดยมีกลไกรับมือกับข้อมูลที่ระบุชื่อไม่ได้ (Unknown / Uncultured) อย่างชาญฉลาด

**โค้ดของฟังก์ชัน:**
```python
    def _parse_taxonomy(self, target_rank):
        """
        จัดการชื่อ Taxonomy ตามเงื่อนไข:
        1. รวมกลุ่มข้อมูลตาม Rank ที่เลือก
        2. ถ้าชื่อเป็น unknown/uncultured ให้ใช้ Rank สูงกว่าแทน
        """
        target_rank = target_rank.lower()
        rank_map = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
        full_rank_map = {'domain': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6}
        
        # ตรวจสอบและแปลง Rank ที่ผู้ใช้พิมพ์เป็น Index ตัวเลข
        if target_rank in full_rank_map:
            rank_idx = full_rank_map[target_rank]
        else:
            rank_idx = rank_map.get(target_rank[0], 5)
            
        taxa_mapping = {}
        for col in self.taxa_cols:
            parts = col.split(';')
            name = 'Unknown'
            # วนลูปถอยหลัง จาก Rank ที่ต้องการ ไปจนถึง Domain
            for i in range(min(rank_idx, len(parts) - 1), -1, -1):
                part = parts[i].strip()
                # ตัด prefix (เช่น g__) ออกมาเพื่อตรวจเช็คความถูกต้องของชื่อ
                val = part[3:] if len(part) > 3 and part[1:3] == '__' else part
                
                # ตรวจสอบว่าเป็นชื่อที่ใช้งานได้หรือไม่
                if val.lower() not in ['', 'incertae_sedis', 'unknown', 'uncultured', 'unassigned'] and not part.endswith('__'):
                    name = part
                    break
            taxa_mapping[col] = name
            
        return taxa_mapping
```

### คำอธิบายการทำงาน (Step-by-Step Logic):

1. **การจับคู่ระดับชั้น (Rank Mapping):**
   * เริ่มต้นด้วยการสร้าง Dictionary (`rank_map` และ `full_rank_map`) เพื่อแปลงชื่อระดับชั้นที่ผู้ใช้ป้อนมา (เช่น `'genus'`, `'g'`, หรือ `'phylum'`) ให้กลายเป็น **ตัวเลข Index** (Domain=0, Phylum=1, ... Species=6) เพื่อนำไปชี้ตำแหน่งในสายอักขระ Taxonomy

2. **การแยกสายอักขระ (String Splitting):**
   * วนลูปรายชื่อคอลัมน์ Taxonomy (เช่น `d__Bacteria;p__Chloroflexota;...`) แล้วใช้ `.split(';')` เพื่อหั่นข้อความออกเป็นส่วนๆ เก็บไว้ในลิสต์ชื่อ `parts` 

3. **กลไกการถอยหลังหาชื่อที่ระบุได้ (Fallback Mechanism):**
   * นี่คือหัวใจสำคัญของฟังก์ชันนี้! โค้ดจะใช้ลูป `for` แบบย้อนกลับ (`step=-1`) โดยเริ่มจาก `rank_idx` (ระดับที่ผู้ใช้ต้องการ) ถอยหลังกลับขึ้นไปยังระดับที่ใหญ่กว่า
   * ในแต่ละรอบ โค้ดจะดึงชื่อมาลบช่องว่าง และแยก Prefix ออก (เช่น ตัด `g__` ออกไปเพื่อเช็คข้อความเพียวๆ ในตัวแปร `val`)
   * **เช็คเงื่อนไขคำต้องห้าม (`if val.lower() not in [...]`):** ระบบจะตรวจสอบว่าชื่อนั้น **ไม่ใช่** คำว่างเปล่า, คำว่า `uncultured`, `unknown`, `incertae_sedis`, `unassigned` หรือคำที่ลงท้ายด้วย `__` 
   * หากชื่อนั้นเป็น "ชื่อจุลินทรีย์ที่แท้จริง" โค้ดจะจดจำค่านั้นและทำการ `break` หยุดค้นหาทันที
   * **กรณีฉุกเฉิน:** หากชื่อในชั้นนั้นดันเป็น `uncultured` หรือหาไม่พบ ลูปจะไม่ถูก `break` และมันจะวนรอบถัดไปเพื่อดึงชื่อใน Rank ที่สูงกว่า (เช่น ดึงชื่อ Family มาใช้แทน Genus) 

4. **การคืนค่าผลลัพธ์ (Return):**
   * ฟังก์ชันจะคืนค่าออกมาเป็นตารางจับคู่ (Dictionary) ระหว่าง "ชื่อเต็มดั้งเดิม" ชี้ไปยัง "ชื่อใหม่ที่ถูกคัดกรองแล้ว" ผลลัพธ์นี้จะถูกส่งต่อให้ Pandas จัดการทำ Group By รวมกลุ่มที่มีชื่อตรงกันเข้าด้วยกันต่อไป

---

## การทำงานของฟังก์ชัน `get_relative_abundance`

ฟังก์ชัน `get_relative_abundance` ทำหน้าที่ **คำนวณค่า Relative Abundance แล้วยุบรวมข้อมูลตามกลุ่ม Metadata** (เช่น Organic / Chemical) เพื่อเตรียมข้อมูลให้พร้อมสำหรับการพล็อต Stacked Bar Plot

**โค้ดของฟังก์ชัน:**
```python
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
```

### คำอธิบายการทำงาน (Step-by-Step Logic):

1. **จับคู่ชื่อ Taxonomy ใหม่ (Taxonomy Mapping):**
   * เรียก `_parse_taxonomy(taxonomic_rank)` เพื่อสร้าง Dictionary ที่จับคู่ "ชื่อ Taxonomy เต็ม" ไปยัง "ชื่อระดับชั้นที่ต้องการ" โดยจัดการ Unknown/Uncultured ให้แล้ว

2. **คำนวณ Relative Abundance ระดับ Sample:**
   * นำตาราง count `df_counts` มาหารด้วยผลรวมของแต่ละแถว (`df_counts.sum(axis=1)`) แล้วคูณ 100 เพื่อให้ได้ค่าเปอร์เซ็นต์
   * ผลลัพธ์คือ แต่ละ Sample (แถว) จะมีผลรวม Relative Abundance ของทุก Taxa เท่ากับ 100%

3. **ยุบรวมคอลัมน์ตาม Taxonomy ใหม่ (Group by Taxa):**
   * ใช้ `.rename(columns=taxa_map)` เปลี่ยนชื่อคอลัมน์เป็นชื่อระดับที่ต้องการ
   * จากนั้น Transpose แล้วทำ `groupby(...).sum()` เพื่อรวมคอลัมน์ที่มีชื่อตรงกัน (เพราะหลายคอลัมน์อาจถูก map เป็นชื่อเดียวกัน)

4. **ยุบรวมตามกลุ่ม Metadata (Group by Farm Type):**
   * เพิ่มคอลัมน์ `'Group'` จาก Metadata คอลัมน์สุดท้าย (Organic/Chemical)
   * ใช้ `.groupby('Group').mean()` เพื่อหาค่าเฉลี่ย Relative Abundance ภายในแต่ละกลุ่ม → ได้ตาราง 1 แถวต่อ 1 กลุ่ม

5. **กรองด้วย Cutoff (Apply Cutoff):**
   * ส่งต่อให้ `_apply_cutoff` กรอง Taxa ที่มีค่า max ต่ำกว่า `percent_cutoff` แล้วยุบรวมเป็น `'Others'`
   * ผลลัพธ์สุดท้ายคือ DataFrame ที่แถวเป็นกลุ่ม (Organic/Chemical) และคอลัมน์เป็น Taxa ที่ผ่านเกณฑ์ พร้อมจัดเรียงตามค่าเฉลี่ยจากมากไปน้อย

---

## การทำงานของฟังก์ชัน `get_abundance_matrix`

ฟังก์ชัน `get_abundance_matrix` ทำหน้าที่ **คำนวณค่า Relative Abundance ในระดับ Sample โดยไม่ยุบรวมตามกลุ่ม Metadata** เพื่อให้ได้เมทริกซ์ที่เหมาะกับการพล็อต Heatmap และ Bubble Plot (ที่แกน X จะเป็นชื่อแต่ละ Sample ไม่ใช่กลุ่ม)

**โค้ดของฟังก์ชัน:**
```python
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
```

### คำอธิบายการทำงาน (Step-by-Step Logic):

1. **จับคู่ชื่อ Taxonomy ใหม่:**
   * เรียก `_parse_taxonomy` เช่นเดียวกับ `get_relative_abundance` เพื่อให้ชื่อ Taxa ถูกจัดให้อยู่ในระดับที่ต้องการ

2. **คำนวณ Relative Abundance ต่อ Sample:**
   * หารด้วยผลรวมแถว แล้วคูณ 100 → แต่ละ Sample รวมกันได้ 100%

3. **ยุบรวมคอลัมน์ที่ชื่อตรงกัน:**
   * Rename + Transpose + `groupby(...).sum()` เพื่อนำคอลัมน์ที่ถูก map เป็นชื่อเดียวกันมารวมกัน

4. **กรองด้วย Cutoff (แบบเก็บราย Sample):**
   * ต่างจาก `_apply_cutoff` ของ Stacked Bar: ที่นี่ **ไม่มีการสร้างกลุ่ม `'Others'`** แต่ใช้ `df_rel.max() >= percent_cutoff` เพื่อเก็บเฉพาะ Taxa ที่มี Sample อย่างน้อยหนึ่งตัวเกินเกณฑ์ที่กำหนด
   * วิธีนี้จะทำให้ Heatmap/Bubble Plot ไม่มีแถว Taxa ที่ค่าต่ำเตี้ยทุก Sample (ช่วยลด noise)

5. **เรียงลำดับ Taxa และคืนค่า:**
   * เรียงคอลัมน์ตามค่าเฉลี่ย Relative Abundance จากมากไปน้อย
   * คืนค่าเป็น DataFrame โดย **rows = Samples** และ **cols = Taxa**

> หมายเหตุ: เปรียบเทียบกับ `get_relative_abundance` — ตัวนั้นทำ `groupby('Group').mean()` แล้ว apply cutoff + เพิ่ม Others ส่วน `get_abundance_matrix` คงระดับ Sample เอาไว้และตัด Taxa ที่ต่ำกว่าเกณฑ์ทิ้งโดยไม่ยุบเป็น Others

---

## การทำงานของฟังก์ชัน `heatmap`

ฟังก์ชัน `heatmap` ทำหน้าที่ **สร้างแผนภาพ Heatmap แสดงรูปแบบความหนาแน่นของจุลินทรีย์ในแต่ละ Sample** โดยรองรับการปรับสเกล (z-score, log) และการจัดกลุ่มแบบ Hierarchical Clustering พร้อมแสดงแถบสีระบุกลุ่ม Metadata ของแต่ละ Sample

**โค้ดของฟังก์ชัน:**
```python
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
        ...
```

### คำอธิบายการทำงาน (Step-by-Step Logic):

1. **เตรียมเมทริกซ์ข้อมูล (Data Preparation):**
   * เรียก `get_abundance_matrix(taxonomic_rank, percent_cutoff)` เพื่อให้ได้ DataFrame ที่ rows = Samples, cols = Taxa
   * เมทริกซ์นี้ผ่านการคำนวณ Relative Abundance, ยุบรวมตามระดับ Taxonomy, และกรอง Taxa ที่ต่ำกว่า Cutoff มาแล้ว

2. **เลือกวิธีปรับสเกล (Scaling Strategy):**
   * **`'log'`** : ใช้ `np.log1p(df)` (คือ `log(1 + x)`) เพื่อบีบช่วงค่าที่กว้างให้อ่านง่าย กันปัญหา log(0) ใช้ colormap `'viridis'`
   * **`'z-score'`** *(ค่าเริ่มต้น)* : คำนวณ z-score ราย Taxon ด้วยสูตร `(x - mean) / std` โดยมีการแทนค่า std ที่เป็น 0 ด้วย 1 เพื่อกัน Division by Zero ใช้ colormap `'vlag'` ที่มีค่ากลางเป็น 0 (เหมาะกับการดู "สูงกว่า/ต่ำกว่าค่าเฉลี่ย")
   * **`None`** : ใช้ค่า Relative Abundance ดิบ ใช้ colormap `'viridis'`

3. **Transpose เพื่อให้ตรงตามสเปคของภาพ:**
   * Heatmap ต้องการ rows = Taxa, cols = Samples จึงทำ `data.T` เพื่อสลับแกน

4. **สร้างแถบสีระบุกลุ่ม Metadata (Column Colors):**
   * ดึงคอลัมน์ Metadata ของแต่ละ Sample (Organic/Chemical) มาจาก `self.df[self.metadata_col]`
   * กำหนด palette เริ่มต้น `Organic = สีเขียว`, `Chemical = สีแดง` หากมีกลุ่มอื่นจะใช้สีจาก `seaborn.color_palette('Set2', ...)` แทน
   * ได้ Series `col_colors` ที่จะนำไปวางเป็นแถบสีเหนือ Heatmap

5. **เลือกโหมดการวาด (Clustering vs ไม่ Clustering):**
   * **`clustering=True`** : ใช้ `sns.clustermap` ซึ่งจะคำนวณ Hierarchical Clustering ทั้งฝั่ง rows (Taxa) และ cols (Samples) แล้วแสดง Dendrogram + จัดเรียง Sample/Taxa ใหม่ตามความใกล้เคียง พร้อมวางแถบ `col_colors` และ Legend ระบุชื่อกลุ่มไว้ด้านบน
   * **`clustering=False`** : ใช้ `sns.heatmap` ธรรมดา จัดเรียงตามลำดับเดิม (Taxa ถูก sort ตาม mean abundance) เหมาะกับการเทียบ Sample ตามลำดับเดิม

6. **ขนาดภาพแบบยืดหยุ่น (Adaptive Figure Size):**
   * ความกว้าง = `max(8, 0.6 * จำนวน Sample + 4)` และความสูง = `max(6, 0.25 * จำนวน Taxa + 2)`
   * ทำให้ภาพไม่อัดแน่นเมื่อมี Sample/Taxa เยอะ และไม่เล็กเกินไปเมื่อมีน้อย

7. **คืนค่า:**
   * ฟังก์ชันคืน object `plt` (pyplot) ผู้ใช้สามารถนำไป `.show()` ดูบนหน้าจอ หรือใช้ `export_heatmap` เพื่อเซฟเป็นไฟล์ภาพในนามสกุลที่ต้องการได้

### ตัวอย่างการใช้งาน:
```python
# Heatmap แบบ z-score พร้อม Hierarchical Clustering (ค่าเริ่มต้น)
plt_fig = visualizer.heatmap(taxonomic_rank='genus', percent_cutoff=1.0,
                             scaling='z-score', clustering=True)
plt_fig.show()

# บันทึก Heatmap แบบ log scale ไม่ทำ Clustering เป็น PNG
visualizer.export_heatmap(taxonomic_rank='phylum', percent_cutoff=2.0,
                          scaling='log', clustering=False,
                          file_format='png', filename='phylum_heatmap')
```
