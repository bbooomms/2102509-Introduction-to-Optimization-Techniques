# Newton Optimization (Rosenbrock Function)

โค้ดนี้เป็นส่วนหนึ่งของการบ้าน/โปรเจกต์วิชา Optimization  
จัดทำและเขียนโดย **Phubet Suwanno (Boom)**  
โดยใช้ Newton’s method ร่วมกับ line search (Golden Section Search) เพื่อหาค่าต่ำสุดของฟังก์ชัน Rosenbrock

---

## ⚠️ หมายเหตุ (Note)
งานนี้จัดทำโดยผมเองทั้งหมด ยังไม่ได้รับการตรวจหรือรีวิวจากอาจารย์หรือผู้สอน อาจมีข้อผิดพลาดทางตรรกะ สมการ หรือการคำนวณบางส่วนได้ หากพบข้อผิดพลาดหรือข้อเสนอแนะ สามารถแจ้งได้เลยครับ ยินดีรับฟังและแก้ไขเสมอ 😊  

---

## 🇹🇭 คำอธิบายโดยย่อ (Thai)
โค้ดนี้ประกอบด้วย:
- `Newton.m` — ฟังก์ชันหลักสำหรับ Newton’s method  
- `Fcn.m` — นิยามของ Rosenbrock function พร้อม gradient และ Hessian  
- `golden.m` — โค้ด Golden Section Search ที่นำมาจากการบ้านครั้งก่อน  
- `run_cases.m` — สคริปต์สำหรับรันการทดลอง 3 เคสที่กำหนด  

พารามิเตอร์ที่ใช้:  
\[
\varepsilon_{\text{rel}} = 0.05, \quad \varepsilon_{\text{abs}} = 10^{-3}
\]

---

## 🇬🇧 English Version
This project was written and implemented by **Phubet Suwanno (Boom)** as part of an optimization assignment. It applies **Newton’s method with line search (Golden Section Search)** to minimize the **Rosenbrock function**.

> ⚠️ This work has not yet been reviewed or verified.  
> There may be minor errors or imperfections in logic, equations, or implementation.  
> Suggestions, corrections, or improvements are warmly welcomed! 😊

**Files included:**
- `Newton.m` — main Newton optimization algorithm  
- `Fcn.m` — Rosenbrock function, gradient, and Hessian  
- `golden.m` — Golden Section Search code from the previous homework  
- `run_cases.m` — test script for three initial points  

Parameters used:
\[
\varepsilon_{\text{rel}} = 0.05, \quad \varepsilon_{\text{abs}} = 10^{-3}
\]

---

## 🧠 Remarks
The code is written for learning purposes and may not be fully optimized for performance. Feel free to modify, reuse, or improve it for your own experiments.

---

**Author:**  
Phubet Suwanno (Boom)  
Faculty of Engineering, Chulalongkorn University  
