# تبدأ منين؟ (النسخة النظيفة — PNG فقط)

1. **النتائج:** `docs/RESULTS_AR.md` للملخص العربي أو `docs/RESULTS_EN.md` للإنجليزي.
2. **الصور:** `figures/main/` فيه 5 صور نتائج PNG، و`figures/supplementary/` فيه 21 صورة Single-cell PNG (تشمل نسخ UMAP الواضحة).
3. **فهرس كامل للرسومات:** `docs/FIGURE_INDEX.md`؛ والأسماء القديمة والجديدة في `docs/FIGURE_FILE_MAPPING.csv`.
4. **GitHub:** المستودع موجود بالفعل على GitHub. **لا تحذفه**؛ استبدل الملفات ببرنامج GitHub Desktop حسب `docs/GITHUB_REPLACE_GUIDE_AR.md`.
5. **مجموعة التحليل:** `primary_07e` معتمدة للتحليل الأساسي فقط؛ القيود في `docs/ANALYSIS_SET_APPROVAL.md` و`docs/METHODS_AND_LIMITATIONS.md`.
6. **اختياري لإعادة الإنتاج:** تصدير exact per-cell masks يتطلب تشغيل `scripts/00_export_and_verify_primary_07e_masks.R` على RStudio مع ملفات المشروع الأصلية.

**مهم:** الارتفاع الوصفي لـ SPP1 وFAP في الورم لا يثبت SPP1–FAP–CD8 coupling على مستوى المرضى.
