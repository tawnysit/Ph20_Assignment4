IMAGES= analytic_solution.jpg explicit_solution.jpg analytic_vs_explicit.jpg explicit_error.jpg truncation_error.jpg explicit_energy.jpg implicit_solution.jpg analytic_vs_implicit.jpg implicit_vs_explicit.jpg implicit_error.jpg implicit_vs_explicit_error.jpg implicit_energy.jpg analytic_phasespace.jpg implicit_vs_explicit_vs_analytic_phasespace.jpg sympletic_phasespace.jpg analytic_vs_sympletic_phasespace.jpg analytic_vs_sympletic_energy.jpg

pdf_writeup.pdf: $(IMAGES)
	pdflatex main.tex $^ > $@

$(IMAGES):
	python euler.py > $@