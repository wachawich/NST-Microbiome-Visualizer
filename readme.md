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
