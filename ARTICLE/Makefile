# Installed publishers and science from texlive to get algorithmic.sty
# and IEEEtran.cls

FIGURES = DR_gp ICY meantime dag density_param_tau_pi trajectory_tau_pi_sum
PDFFIGS = $(addsuffix .pdf, ${FIGURES})

TARGET = bayesianbaselineestimation
${TARGET}.pdf : ${TARGET}.tex ${PDFFIGS}
	pdflatex ${TARGET}

## variables     : Print selected variables.
.PHONY : variables
variables:
	@echo FIGURES : ${FIGURES}
	@echo PDFFIGS : ${PDFFIGS}

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
