# دليل الرسومات — النسخة المبسّطة PNG فقط

- النتائج الرئيسية: `figures/main/` — خمس صور PNG من جداول R المشتقة، مرسومة بواسطة Python.
- صور الـ Single-cell: `figures/supplementary/` — 21 صورة PNG مأخوذة من تحليل R/Seurat أو من صور QC الأصلية.
- الصور الباهتة القديمة `S02` و`S03` تم استبدالها بنسخ عالية الوضوح من R، بالإضافة لنسخة Tumor/Normal مقسمة.
- أسماء الرسومات ووصفها في `docs/FIGURE_INDEX.md`، ومطابقة القديم والجديد في `docs/FIGURE_FILE_MAPPING.csv`.
- لا توجد ملفات SVG في النسخة دي، ولا توجد صور workflow من الذكاء الاصطناعي.
- بعض الـ clusters في صور Validation ما زالت مسماة على مستوى subclustering المبدئي؛ الـ masks التحليلية معتمدة لكن ده مش إثبات تجريبي للـ cell identities.
- السكريبت `scripts/01_reproduce_figures.py` بيعيد توليد الـ 5 PNG الرئيسية؛ وسكريبت R للرسومات الفرعية ضمن `scripts/02_generate_singlecell_QC_annotation_figures.R`.

## رفع المشروع بدل النسخة القديمة
لا ترفع ZIP على GitHub كملف؛ استبدل الملفات عن طريق GitHub Desktop (Clone → copy new contents → Commit → Push).
الطريقة دي هتمسح الملفات القديمة اللي مش موجودة في النسخة الجديدة وتحافظ على نفس رابط الريبو.
