#!/usr/bin/env python3
"""
Generate PDF report from Markdown using FPDF2.

This script converts VLE_Consulting_Report.md to a professional PDF format.
"""

import os
import re
from fpdf import FPDF

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DOCS_DIR = os.path.join(SCRIPT_DIR, "..", "docs")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
PDF_FILE = os.path.join(DOCS_DIR, "VLE_Consulting_Report.pdf")


class PDFReport(FPDF):
    """Custom PDF class with header and footer."""
    
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=25)
    
    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 9)
            self.set_text_color(128)
            self.cell(0, 10, 'VLE Consulting Report - Ethanol-Water System', align='C')
            self.ln(5)
            self.set_draw_color(200)
            self.line(10, 18, 200, 18)
            self.ln(10)
    
    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(128)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', align='C')
    
    def chapter_title(self, title, level=1):
        if level == 1:
            self.set_font('Helvetica', 'B', 16)
            self.set_text_color(26, 26, 26)
            self.ln(5)
            self.cell(0, 10, title, ln=True)
            self.set_draw_color(51)
            self.set_line_width(0.5)
            self.line(10, self.get_y(), 200, self.get_y())
            self.ln(5)
        elif level == 2:
            self.set_font('Helvetica', 'B', 14)
            self.set_text_color(44, 62, 80)
            self.ln(3)
            self.cell(0, 8, title, ln=True)
            self.set_draw_color(153)
            self.set_line_width(0.3)
            self.line(10, self.get_y(), 150, self.get_y())
            self.ln(3)
        else:
            self.set_font('Helvetica', 'B', 12)
            self.set_text_color(52, 73, 94)
            self.ln(2)
            self.cell(0, 6, title, ln=True)
            self.ln(2)
    
    def body_text(self, text):
        self.set_font('Helvetica', '', 10)
        self.set_text_color(51)
        self.multi_cell(0, 5, text)
        self.ln(2)
    
    def add_image_with_caption(self, img_path, caption):
        if os.path.exists(img_path):
            self.image(img_path, x=25, w=160)
            self.set_font('Helvetica', 'I', 9)
            self.set_text_color(100)
            self.cell(0, 5, caption, ln=True, align='C')
            self.ln(5)


def main():
    print("=" * 50)
    print("PDF Report Generator (FPDF2)")
    print("=" * 50)
    
    pdf = PDFReport()
    pdf.alias_nb_pages()
    pdf.add_page()
    
    # ================== TITLE PAGE ==================
    pdf.ln(60)
    pdf.set_font('Helvetica', 'B', 28)
    pdf.set_text_color(26, 26, 26)
    pdf.cell(0, 15, 'Vapor-Liquid Equilibrium', ln=True, align='C')
    pdf.cell(0, 15, 'Consulting Report', ln=True, align='C')
    
    pdf.ln(10)
    pdf.set_font('Helvetica', '', 16)
    pdf.set_text_color(85)
    pdf.cell(0, 10, 'Non-Ideal Mixing System: Ethanol-Water at 1 atm', ln=True, align='C')
    
    pdf.ln(30)
    pdf.set_font('Helvetica', '', 11)
    pdf.cell(0, 6, 'Prepared for: Management Team', ln=True, align='C')
    pdf.cell(0, 6, 'Fine Chemicals & Biotechnology Materials Co.', ln=True, align='C')
    
    pdf.ln(10)
    pdf.cell(0, 6, 'Prepared by: VLE Consulting Team', ln=True, align='C')
    pdf.cell(0, 6, 'January 2026', ln=True, align='C')
    
    pdf.ln(20)
    pdf.set_font('Helvetica', 'I', 9)
    pdf.set_text_color(128)
    pdf.cell(0, 5, 'Reference: Sandler, S.I., Chemical, Biochemical, and Engineering Thermodynamics', ln=True, align='C')
    pdf.cell(0, 5, 'Chapters 6-10', ln=True, align='C')
    
    # ================== PART I ==================
    pdf.add_page()
    pdf.chapter_title('Part I: Executive Summary', 1)
    
    pdf.chapter_title('The Core Thermodynamic Problem', 2)
    pdf.body_text(
        "The client's distillation and solvent recovery processes have experienced significant "
        "operational failures due to the application of ideal solution assumptions (Raoult's Law) "
        "to the Ethanol-Water binary system. This fundamental modeling error has resulted in:\n\n"
        "  - Energy consumption deviations exceeding 15-20% from design specifications\n"
        "  - Unstable operating temperatures in the distillation column\n"
        "  - Failure to achieve expected phase separation behavior"
    )
    
    pdf.chapter_title('Why Ideal Solution Assumptions Failed', 2)
    pdf.body_text(
        "The Ethanol-Water system exhibits strong positive deviations from ideal behavior due to:\n\n"
        "  1. Hydrogen bonding disruption when ethanol and water molecules mix\n"
        "  2. Molecular size asymmetry between ethanol (C2H5OH) and water (H2O)\n"
        "  3. Self-association differences between the two polar molecules\n\n"
        "These non-ideal interactions result in activity coefficients significantly greater than "
        "unity (gamma >> 1), which Raoult's Law fundamentally cannot capture."
    )
    
    pdf.chapter_title('Primary Engineering Recommendation', 2)
    pdf.body_text(
        "We strongly recommend redesigning the distillation process to incorporate:\n\n"
        "  1. A pressure-swing or extractive distillation scheme to break the azeotrope\n"
        "  2. Recalculated heat exchanger duties using NRTL-based heat of mixing data\n"
        "  3. Operating temperature targets adjusted based on accurate bubble point predictions"
    )
    
    # ================== PART II ==================
    pdf.add_page()
    pdf.chapter_title('Part II: Problem Definition and Theoretical Background', 1)
    
    pdf.chapter_title('Unit Operation Context', 2)
    pdf.body_text(
        "The unit operation under evaluation is a continuous distillation column for solvent "
        "recovery, where:\n\n"
        "  - Feed: Ethanol-Water mixture from upstream fermentation/extraction\n"
        "  - Products: Purified ethanol (overhead) and water (bottoms)\n"
        "  - Operating Pressure: 1 atm (101.3 kPa)"
    )
    
    pdf.chapter_title('Phase Equilibrium Considerations', 2)
    pdf.body_text(
        "This process involves Vapor-Liquid Equilibrium (VLE) where both liquid and vapor "
        "phases coexist at bubble/dew point conditions.\n\n"
        "In non-ideal liquid mixtures, the fugacity of component i in the liquid phase is:\n"
        "   fi(L) = xi * gamma_i * fi(0,L)\n\n"
        "For ideal solutions, gamma_i = 1 (Raoult's Law). However, for Ethanol-Water, "
        "activity coefficients range from 1.5 to 4.5."
    )
    
    # ================== PART III ==================
    pdf.add_page()
    pdf.chapter_title('Part III: NRTL Model and Phase Equilibrium Analysis', 1)
    
    pdf.chapter_title('NRTL Model Framework', 2)
    pdf.body_text(
        "The Non-Random Two-Liquid (NRTL) model, developed by Renon and Prausnitz (1968), "
        "is based on the Local Composition concept. Key parameters:\n\n"
        "  - tau_ij: Energy parameter representing interaction energy differences\n"
        "  - alpha_ij: Non-randomness factor (typically 0.2-0.47)\n\n"
        "For Ethanol-Water: alpha_12 = alpha_21 = 0.3, tau_12 = 0.88, tau_21 = 1.45"
    )
    
    # Add T-x-y diagram
    txy_path = os.path.join(FIGURES_DIR, "TxY_Ethanol_Water_1atm.png")
    if os.path.exists(txy_path):
        pdf.ln(5)
        pdf.add_image_with_caption(txy_path, "Figure 1: T-x-y Diagram comparing Ideal, Van Laar, and NRTL models")
    
    pdf.add_page()
    # Add y-x diagram
    yx_path = os.path.join(FIGURES_DIR, "Yx_Ethanol_Water_1atm.png")
    if os.path.exists(yx_path):
        pdf.add_image_with_caption(yx_path, "Figure 2: y-x Diagram showing azeotropic behavior")
    
    pdf.chapter_title('Key Findings', 2)
    pdf.body_text(
        "  - Maximum deviation at x1 = 0.3-0.5 (ethanol-lean compositions)\n"
        "  - Temperature prediction error: up to 5-8 deg C in mid-composition range\n"
        "  - Azeotrope detected at x1 = 0.89 (95.6 wt% ethanol)"
    )
    
    # ================== PART IV ==================
    pdf.add_page()
    pdf.chapter_title('Part IV: Excess Properties and Energy Analysis', 1)
    
    # Add G^E diagram
    ge_path = os.path.join(FIGURES_DIR, "GE_vs_x_Ethanol_Water_1atm.png")
    if os.path.exists(ge_path):
        pdf.add_image_with_caption(ge_path, "Figure 3: Excess Gibbs Energy vs Composition")
    
    # Add H^E diagram
    he_path = os.path.join(FIGURES_DIR, "Hmix_vs_x_Ethanol_Water_1atm.png")
    if os.path.exists(he_path):
        pdf.add_image_with_caption(he_path, "Figure 4: Excess Enthalpy (Heat of Mixing) vs Composition")
    
    pdf.chapter_title('Thermal Behavior Analysis', 2)
    pdf.body_text(
        "  - G^E sign: Positive throughout - confirms positive deviation from Raoult's Law\n"
        "  - Maximum G^E approximately 0.8-1.0 kJ/mol at x1 = 0.4\n"
        "  - Heat effects must be incorporated into heat exchanger sizing"
    )
    
    # ================== PART V ==================
    pdf.add_page()
    pdf.chapter_title('Part V: Consulting Recommendations', 1)
    
    pdf.chapter_title('Process Feasibility Assessment', 2)
    pdf.body_text(
        "Current operating conditions are NOT feasible for achieving high-purity ethanol "
        "(>96 mol%) through simple distillation due to the azeotrope at x1 = 0.89."
    )
    
    pdf.chapter_title('Engineering Recommendations', 2)
    pdf.body_text(
        "1. PRESSURE-SWING DISTILLATION\n"
        "   Operate two columns at different pressures (1 atm and 0.2 atm)\n"
        "   Estimated additional capital: 30-40% increase\n\n"
        "2. EXTRACTIVE DISTILLATION\n"
        "   Add a high-boiling entrainer (e.g., ethylene glycol)\n"
        "   Recommended for large-scale, continuous production\n\n"
        "3. HYBRID PERVAPORATION\n"
        "   Distillation to ~90 mol%, then membrane pervaporation\n"
        "   Lowest energy consumption option for fuel-grade ethanol"
    )
    
    pdf.ln(10)
    pdf.set_font('Helvetica', 'B', 11)
    pdf.set_text_color(26, 26, 26)
    pdf.cell(0, 6, 'VERDICT: Process redesign required using NRTL-based VLE predictions', ln=True, align='C')
    
    # Save PDF
    pdf.output(PDF_FILE)
    
    print(f"\nPDF generated successfully!")
    print(f"Location: {PDF_FILE}")
    print("=" * 50)


if __name__ == "__main__":
    main()
