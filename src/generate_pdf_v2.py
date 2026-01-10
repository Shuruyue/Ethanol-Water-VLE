#!/usr/bin/env python3
"""
Generate PDF report V2 from Markdown using FPDF2.
Target: <= 20 pages with complete mathematical models.
"""

import os
from fpdf import FPDF

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCS_DIR = os.path.join(SCRIPT_DIR, "..", "docs")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
PDF_FILE = os.path.join(DOCS_DIR, "VLE_Consulting_Report_V2.pdf")


class PDFReportV2(FPDF):
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=20)
    
    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(100)
            self.cell(95, 8, 'VLE Consulting Report V2', align='L')
            self.cell(95, 8, f'Page {self.page_no()}', align='R')
            self.ln(3)
            self.set_draw_color(180)
            self.line(10, 15, 200, 15)
            self.ln(8)
    
    def footer(self):
        self.set_y(-12)
        self.set_font('Helvetica', 'I', 7)
        self.set_text_color(128)
        self.cell(0, 5, 'ChemSep Validated | thermo Library | January 2026', align='C')
    
    def chapter_title(self, title, level=1):
        if level == 1:
            self.set_font('Helvetica', 'B', 14)
            self.set_text_color(20, 60, 100)
            self.ln(4)
            self.cell(0, 8, title, ln=True)
            self.set_draw_color(20, 60, 100)
            self.set_line_width(0.6)
            self.line(10, self.get_y(), 200, self.get_y())
            self.ln(4)
        elif level == 2:
            self.set_font('Helvetica', 'B', 11)
            self.set_text_color(40, 80, 120)
            self.ln(2)
            self.cell(0, 6, title, ln=True)
            self.ln(1)
        else:
            self.set_font('Helvetica', 'B', 10)
            self.set_text_color(60, 60, 60)
            self.cell(0, 5, title, ln=True)
    
    def body_text(self, text):
        self.set_font('Helvetica', '', 9)
        self.set_text_color(40)
        self.multi_cell(0, 4.5, text)
        self.ln(1)
    
    def equation(self, eq_text, label=None):
        self.set_fill_color(245, 245, 250)
        self.set_font('Courier', '', 9)
        self.set_text_color(30)
        self.ln(1)
        self.multi_cell(0, 5, eq_text, fill=True, align='C')
        if label:
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(100)
            self.cell(0, 4, label, ln=True, align='R')
        self.ln(1)
    
    def add_table(self, headers, data, col_widths=None):
        if col_widths is None:
            col_widths = [190 / len(headers)] * len(headers)
        self.set_font('Helvetica', 'B', 8)
        self.set_fill_color(20, 60, 100)
        self.set_text_color(255)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 6, h, border=1, fill=True, align='C')
        self.ln()
        self.set_font('Helvetica', '', 8)
        self.set_text_color(40)
        fill = False
        for row in data:
            self.set_fill_color(245, 248, 250) if fill else self.set_fill_color(255)
            for i, c in enumerate(row):
                self.cell(col_widths[i], 5, str(c), border=1, fill=True, align='C')
            self.ln()
            fill = not fill
        self.ln(2)
    
    def add_image(self, img_path, caption, width=140):
        if os.path.exists(img_path):
            self.image(img_path, x=(210-width)/2, w=width)
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(80)
            self.cell(0, 4, caption, ln=True, align='C')
            self.ln(3)


def main():
    print("=" * 60)
    print("VLE Consulting Report V2 - PDF Generator")
    print("=" * 60)
    
    pdf = PDFReportV2()
    pdf.alias_nb_pages()
    
    # TITLE PAGE
    pdf.add_page()
    pdf.ln(40)
    pdf.set_font('Helvetica', 'B', 24)
    pdf.set_text_color(20, 60, 100)
    pdf.cell(0, 12, 'Vapor-Liquid Equilibrium', ln=True, align='C')
    pdf.cell(0, 12, 'Consulting Report', ln=True, align='C')
    pdf.ln(5)
    pdf.set_font('Helvetica', '', 14)
    pdf.set_text_color(80)
    pdf.cell(0, 8, 'Non-Ideal Binary System: Ethanol-Water at 1 atm', ln=True, align='C')
    pdf.ln(10)
    pdf.set_font('Helvetica', 'B', 18)
    pdf.set_text_color(200, 50, 50)
    pdf.cell(0, 10, 'VERSION 2.0', ln=True, align='C')
    pdf.ln(15)
    pdf.set_font('Helvetica', '', 10)
    pdf.set_text_color(60)
    pdf.cell(0, 5, 'Prepared for: Management Team', ln=True, align='C')
    pdf.cell(0, 5, 'Fine Chemicals & Biotechnology Materials Co.', ln=True, align='C')
    pdf.ln(5)
    pdf.cell(0, 5, 'Prepared by: VLE Consulting Team | January 2026', ln=True, align='C')
    pdf.ln(10)
    pdf.set_font('Helvetica', 'I', 8)
    pdf.set_text_color(100)
    pdf.cell(0, 4, 'Reference: Sandler, Chemical Engineering Thermodynamics, 5th Ed.', ln=True, align='C')
    pdf.cell(0, 4, 'Parameters: ChemSep Database (Validated)', ln=True, align='C')
    
    # TABLE OF CONTENTS
    pdf.add_page()
    pdf.chapter_title('Table of Contents', 1)
    toc = [("1. Executive Summary", 3), ("2. Thermodynamic Fundamentals", 4),
           ("3. Mathematical Framework", 5), ("4. NRTL Model", 7),
           ("5. Phase Equilibrium Analysis", 9), ("6. Excess Properties", 12),
           ("7. Engineering Recommendations", 14), ("Appendix", 16)]
    pdf.set_font('Helvetica', '', 10)
    for t, p in toc:
        pdf.set_text_color(40)
        pdf.cell(160, 6, t)
        pdf.set_text_color(100)
        pdf.cell(20, 6, str(p), align='R')
        pdf.ln()
    
    # 1. EXECUTIVE SUMMARY
    pdf.add_page()
    pdf.chapter_title('1. Executive Summary', 1)
    pdf.chapter_title('Problem Statement', 2)
    pdf.body_text("The client's Ethanol-Water distillation has failed due to ideal solution assumptions. "
                  "This report uses validated NRTL activity coefficient model for accurate predictions.")
    pdf.chapter_title('Key Findings', 2)
    pdf.add_table(['Parameter', 'Ideal', 'NRTL', 'Impact'],
        [['Azeotrope', 'None', '78.1C, x1=0.874', 'Design limit'],
         ['Max gamma_1', '1.0', '5.7', '470% error'],
         ['Bubble T error', '-', 'Up to 8C', 'Control fail'],
         ['Max G^E', '0', '925 J/mol', 'Non-ideal'],
         ['Max H^E', '0', '663 J/mol', 'Heat effects']], [45, 35, 55, 55])
    pdf.chapter_title('Recommendation', 2)
    pdf.body_text("REDESIGN REQUIRED: Implement pressure-swing or extractive distillation "
                  "to break azeotrope at x1=0.874 (95.6 wt% ethanol).")
    
    # 2. THERMODYNAMIC FUNDAMENTALS
    pdf.add_page()
    pdf.chapter_title('2. Thermodynamic Fundamentals', 1)
    pdf.chapter_title('2.1 Phase Equilibrium Criterion', 2)
    pdf.body_text("At equilibrium between liquid (L) and vapor (V) phases:")
    pdf.equation("f_i^(L) = f_i^(V)    for all components i", "(Eq. 2.1)")
    pdf.chapter_title('2.2 Fugacity Expressions', 2)
    pdf.body_text("Liquid phase with activity coefficient:")
    pdf.equation("f_i^(L) = x_i * gamma_i * f_i^(0,L)", "(Eq. 2.2)")
    pdf.body_text("Vapor phase (ideal gas, low P):")
    pdf.equation("f_i^(V) = y_i * P", "(Eq. 2.3)")
    pdf.chapter_title("2.3 Modified Raoult's Law", 2)
    pdf.equation("y_i * P = x_i * gamma_i * P_i^sat(T)", "(Eq. 2.4)")
    pdf.body_text("For Ethanol-Water, gamma_1 = 1.0 to 5.7 due to H-bond disruption.")
    
    # 3. MATHEMATICAL FRAMEWORK
    pdf.add_page()
    pdf.chapter_title('3. Mathematical Framework', 1)
    pdf.chapter_title('3.1 Excess Gibbs Energy', 2)
    pdf.equation("G^E = RT * sum(x_i * ln(gamma_i))", "(Eq. 3.1)")
    pdf.body_text("Positive G^E -> Positive deviation (gamma > 1) -> Azeotrope possible")
    pdf.chapter_title('3.2 Activity Coefficient', 2)
    pdf.equation("ln(gamma_i) = (1/RT) * [d(nG^E)/dn_i]", "(Eq. 3.2)")
    pdf.chapter_title('3.3 Excess Enthalpy', 2)
    pdf.equation("H^E = -R*T^2 * sum(x_i * d(ln gamma_i)/dT)", "(Eq. 3.3)")
    pdf.body_text("Temperature-dependent parameters essential for H^E calculation.")
    pdf.chapter_title('3.4 Antoine Equation', 2)
    pdf.equation("log10(P^sat) = A - B/(T + C)  [T in C, P in mmHg]", "(Eq. 3.4)")
    pdf.add_table(['Component', 'A', 'B', 'C', 'Source'],
        [['Ethanol', '8.20417', '1642.89', '230.3', 'NIST'],
         ['Water', '8.07131', '1730.63', '233.426', 'NIST']], [35, 35, 35, 30, 55])
    
    # 4. NRTL MODEL
    pdf.add_page()
    pdf.chapter_title('4. NRTL Model Implementation', 1)
    pdf.chapter_title('4.1 Model Origin', 2)
    pdf.body_text("NRTL (Renon & Prausnitz, 1968) uses local composition concept: "
                  "molecular distribution differs from bulk due to interaction energies.")
    pdf.chapter_title('4.2 NRTL Equations', 2)
    pdf.equation("G^E/RT = x1*x2 * [tau21*G21/(x1+x2*G21) + tau12*G12/(x2+x1*G12)]", "(Eq. 4.1)")
    pdf.equation("ln(gamma_1) = x2^2 * [tau21*(G21/(x1+x2*G21))^2 + tau12*G12/(x2+x1*G12)^2]", "(Eq. 4.2)")
    pdf.equation("G_ij = exp(-alpha_ij * tau_ij)    tau_ij = b_ij/T", "(Eq. 4.3)")
    pdf.chapter_title('4.3 Parameter Interpretation', 2)
    pdf.add_table(['Parameter', 'Symbol', 'Meaning'],
        [['Energy', 'tau_ij', 'Interaction energy difference'],
         ['Non-randomness', 'alpha_ij', 'Molecular distribution (0.2-0.47)']], [40, 35, 115])
    pdf.chapter_title('4.4 ChemSep Validated Parameters', 2)
    pdf.body_text("Source: thermo/Interaction Parameters/ChemSep/nrtl.json\n"
                  "Validation: tests/test_nrtl.py::test_NRTL_chemsep")
    pdf.add_table(['Parameter', 'Value', 'Description'],
        [['alpha_12=alpha_21', '0.2937', 'Symmetric'],
         ['b_12', '-29.167 K', 'Ethanol->Water'],
         ['b_21', '+624.868 K', 'Water->Ethanol']], [50, 45, 95])
    pdf.body_text("At T=343.15K: tau_12=-0.085, tau_21=+1.821")
    
    # 5. PHASE EQUILIBRIUM
    pdf.add_page()
    pdf.chapter_title('5. Phase Equilibrium Analysis', 1)
    pdf.chapter_title('5.1 Bubble Point Algorithm', 2)
    pdf.body_text("1. Guess T\n2. Calculate P_sat from Antoine\n3. Calculate tau, G from NRTL\n"
                  "4. Calculate gamma from NRTL\n5. P_calc = x1*gamma1*P1_sat + x2*gamma2*P2_sat\n"
                  "6. Iterate until |P_calc - P| < eps\n7. y1 = x1*gamma1*P1_sat/P")
    pdf.chapter_title('5.2 T-x-y Diagram', 2)
    pdf.add_image(os.path.join(FIGURES_DIR, "TxY_Ethanol_Water_1atm.png"),
                  "Figure 1: T-x-y at 1 atm - Ideal vs NRTL")
    pdf.body_text("NRTL predicts minimum-boiling azeotrope; Ideal shows no azeotrope.")
    
    pdf.add_page()
    pdf.chapter_title('5.3 y-x Diagram', 2)
    pdf.add_image(os.path.join(FIGURES_DIR, "Yx_Ethanol_Water_1atm.png"),
                  "Figure 2: y-x showing azeotrope")
    pdf.chapter_title('5.4 Azeotrope', 2)
    pdf.add_table(['Property', 'Value', 'Implication'],
        [['Type', 'Minimum-boiling', 'T < pure components'],
         ['Temperature', '78.1C', 'Below ethanol NBP'],
         ['Composition', 'x1=0.874', '95.6 wt% ethanol'],
         ['Rel. volatility', 'alpha=1.0', 'No simple distillation']], [45, 50, 95])
    
    # 6. EXCESS PROPERTIES
    pdf.add_page()
    pdf.chapter_title('6. Excess Property Analysis', 1)
    pdf.chapter_title('6.1 Excess Gibbs Energy', 2)
    pdf.add_image(os.path.join(FIGURES_DIR, "GE_vs_x_Ethanol_Water_1atm.png"),
                  "Figure 3: G^E vs composition")
    pdf.add_table(['Property', 'Value', 'Composition'],
        [['Max G^E', '924.55 J/mol', 'x1=0.425'],
         ['G^E at azeotrope', '~500 J/mol', 'x1=0.874']], [60, 60, 70])
    pdf.chapter_title('6.2 Excess Enthalpy', 2)
    pdf.add_image(os.path.join(FIGURES_DIR, "Hmix_vs_x_Ethanol_Water_1atm.png"),
                  "Figure 4: H^E vs composition")
    pdf.add_table(['Property', 'Value', 'Note'],
        [['Max H^E', '663.22 J/mol', 'x1=0.342'],
         ['Behavior', 'Endothermic', 'At bubble T']], [50, 55, 85])
    pdf.chapter_title('6.3 Engineering Impact', 2)
    pdf.add_table(['Effect', 'Magnitude', 'Impact'],
        [['Heat of mixing', '+663 J/mol', 'Heat input needed'],
         ['Non-ideal VLE', 'gamma up to 5.7', 'Larger reflux'],
         ['Azeotrope', 'x1=0.874', 'Max 95.6 wt%']], [50, 50, 90])
    
    # 7. RECOMMENDATIONS
    pdf.add_page()
    pdf.chapter_title('7. Engineering Recommendations', 1)
    pdf.chapter_title('7.1 Azeotrope-Breaking Options', 2)
    pdf.chapter_title('A: Pressure-Swing Distillation', 3)
    pdf.add_table(['Pressure', 'Azeotrope x1', 'T'],
        [['1.0 atm', '0.874', '78.1C'],
         ['0.2 atm', '~0.91', '~52C']], [63, 63, 64])
    pdf.body_text("Two columns at different P can cross azeotrope.")
    pdf.chapter_title('B: Extractive Distillation', 3)
    pdf.body_text("Add entrainer (e.g., ethylene glycol) to break azeotrope.")
    pdf.chapter_title('C: Pervaporation Hybrid', 3)
    pdf.body_text("Distill to 90%, then membrane pervaporation for dehydration.")
    pdf.chapter_title('7.2 Action Summary', 2)
    pdf.add_table(['Priority', 'Action', 'Outcome'],
        [['1', 'Update to NRTL model', 'Accurate predictions'],
         ['2', 'Redesign for azeotrope', '>95.6 wt% ethanol'],
         ['3', 'Recalculate heat duties', 'Proper sizing']], [30, 80, 80])
    pdf.ln(5)
    pdf.set_font('Helvetica', 'B', 11)
    pdf.set_text_color(200, 50, 50)
    pdf.cell(0, 8, 'VERDICT: Process redesign required', ln=True, align='C')
    
    # APPENDIX
    pdf.add_page()
    pdf.chapter_title('Appendix: Technical Data', 1)
    pdf.chapter_title('A.1 Pure Component Properties', 2)
    pdf.add_table(['Property', 'Ethanol', 'Water'],
        [['CAS', '64-17-5', '7732-18-5'],
         ['MW', '46.07 g/mol', '18.02 g/mol'],
         ['NBP', '78.4C', '100.0C'],
         ['Tc', '513.9 K', '647.1 K']], [60, 65, 65])
    pdf.chapter_title('A.2 NRTL Parameters', 2)
    pdf.add_table(['Parameter', 'Value', 'Definition'],
        [['alpha', '0.2937', 'Non-randomness'],
         ['a_12, a_21', '0.0', 'No constant term'],
         ['b_12', '-29.167 K', 'Ethanol->Water'],
         ['b_21', '+624.868 K', 'Water->Ethanol']], [45, 55, 90])
    pdf.chapter_title('A.3 Validation', 2)
    pdf.add_table(['Test', 'gamma_1', 'gamma_2'],
        [['T=343K, x1=0.252', '1.985', '1.146'],
         ['T=343K, x1=0.0', '5.66', '1.00'],
         ['T=343K, x1=1.0', '1.00', '2.67']], [70, 60, 60])
    pdf.chapter_title('A.4 References', 2)
    pdf.set_font('Helvetica', '', 8)
    for r in ["1. Sandler (2017) Chemical Engineering Thermodynamics, 5th Ed.",
              "2. Renon & Prausnitz (1968) AIChE J. 14(1):135-144",
              "3. ChemSep nrtl.json (thermo library)",
              "4. NIST Chemistry WebBook",
              "5. thermo Python library (github.com/CalebBell/thermo)"]:
        pdf.cell(0, 4, r, ln=True)
    
    # SAVE
    pdf.output(PDF_FILE)
    n = pdf.page_no()
    print(f"\nPDF V2 generated: {n} pages")
    print(f"Target: <=20 | Status: {'PASS' if n<=20 else 'FAIL'}")
    print(f"Location: {PDF_FILE}")
    print("=" * 60)


if __name__ == "__main__":
    main()
