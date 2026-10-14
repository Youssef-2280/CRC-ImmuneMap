# استبدال ملفات GitHub بدون حذف الـ Repository

**ملحوظة:** المستودع الحالي Private واسمُه `CRC-ImmuneMap`. المطلوب استبدال محتواه، مش حذف المستودع نفسه.

1. نزّل **GitHub Desktop** من https://desktop.github.com/ وثبّته وسجّل الدخول بنفس حساب GitHub.
2. من GitHub Desktop اختار **File → Clone repository**، ثم اختار `Youssef-2280/CRC-ImmuneMap` واضغط **Clone**.
3. من GitHub Desktop اختار **Repository → Show in Explorer**. هيفتح مجلد الريبو المحلي.
4. احتفظ بمجلد `.git` المخفي الموجود جوه مجلد الريبو، **ومتلمسوش**. احذف ملفات ومجلدات النسخة القديمة ما عدا `.git`.
5. فك ضغط `CRC_ImmuneMap_CLEAN_PNG_v5.zip`، وادخل جوه المجلد اللي فيه `README.md` و`figures/`، وانسخ **محتوياته** إلى مجلد الريبو المحلي.
6. ارجع GitHub Desktop → تبويب **Changes**. لازم تشوف الملفات القديمة اتحذفت (Deleted)، الجديدة اتضافت (Added)، وبعض الملفات اتعدلت (Modified). افحص إن مفيش `.rds` كبيرة أو raw data.
7. اكتب في Summary: `Replace old figures with clear PNG-only release`، واضغط **Commit to main** ثم **Push origin**.
8. افتح GitHub من المتصفح: المفروض تشوف `figures/main/` و`figures/supplementary/` فقط، وملف `docs/FIGURE_INDEX.md`.

لو مش مرتاح تمسح ملفات يدويًا، متعملش الخطوة 4 من نفسك؛ اطلب مني أمشي معاك ضغطة بضغطة في GitHub Desktop.

**تحذير:** حذف كل المستودع من GitHub Settings → Danger Zone مختلف تمامًا ومش مطلوب.
