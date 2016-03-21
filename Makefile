# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R
RFRESH=R --no-init-file --no-environ --no-restore -q
# --vanilla doesn't work with R_PROFILE (because of --no-site-info)

interactive:
	@$(LOAD) R -q --no-save

interactive-emacs:
	@$(LOAD) emacs -nw -f R

.devtools:
	@$(LOAD) $(RFRESH) -e "library(devtools); devtools:::$(FUNC)()"

dependencies: FUNC=install_deps
test: FUNC=test
check: FUNC=check
document: FUNC=document
vignette: FUNC=build_vignettes
clean: FUNC=clean_vignettes
#build: FUNC=build # Can be re-enabled as needed
#build: .devtools
dependencies test check document vignette clean: .devtools
