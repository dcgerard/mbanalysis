nc = 14
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

## Null simulation plots
null_sim_figs = ./output/sims/plots/qq_chisq_g.pdf \
                ./output/sims/plots/qq_chisq_gl.pdf \
                ./output/sims/plots/qq_lrt_g.pdf \
                ./output/sims/plots/qq_lrt_gl.pdf \
                ./output/sims/plots/qq_polymapr_g.pdf \
                ./output/sims/plots/qq_polymapr_gl.pdf \
                ./output/sims/plots/box_lbf_g.pdf \
                ./output/sims/plots/box_lbf_gl.pdf

## Alt sims figs
alt_sim_figs = ./output/sims/plots/t1e.tex \
               ./output/sims/plots/alt_p_box.pdf \
               ./output/sims/plots/alt_lbf_box.pdf

## Blueberry figs
blue_figs = ./output/blue/plots/pairs.pdf \
            ./output/blue/plots/polystrong.pdf \
            ./output/blue/plots/lrtstrong.pdf \
            ./output/blue/plots/tensnps.tex

## Prior sensitivity figs
bayes_figs = ./output/sims/plots/box_lbf_p_g_20.pdf \
                      ./output/sims/plots/box_lbf_p_g_200.pdf \
                      ./output/sims/plots/box_lbf_p_gl_20.pdf \
                      ./output/sims/plots/box_lbf_p_gl_200.pdf \
                      ./output/sims/plots/alt_lbf_box_p.pdf

.PHONY : all
all : sims blue hyp

## Simulations ----

.PHONY : sims
sims : $(null_sim_figs) $(alt_sim_figs) ./output/sims/plots/alpha_ests.pdf $(bayes_figs)

$(null_sim_figs) : ./analysis/null_plots.R ./output/sims/gsims.csv ./output/sims/glsims.csv
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

$(alt_sim_figs) : ./analysis/alt_plots.R ./output/sims/g_altsims.csv ./output/sims/gl_altsims.csv
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

./output/sims/gsims.csv : ./analysis/g_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/sims/g_altsims.csv : ./analysis/g_sims_alt.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/sims/gl_altsims.csv : ./analysis/gl_sims_alt.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/sims/glsims.csv :  ./analysis/gl_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/sims/plots/alpha_ests.pdf : ./analysis/null_plot_ests.R ./output/sims/gsims.csv ./output/sims/glsims.csv
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

$(bayes_figs) : ./analysis/prior_sensitivity.R ./output/sims/gsims.csv ./output/sims/glsims.csv ./output/sims/g_altsims.csv ./output/sims/gl_altsims.csv
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

## Blueberries ----

.PHONY : blue
blue : $(blue_figs)

./output/blue/bluefits.RDS : ./analysis/blue_up.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/blue/blue_df.csv : ./analysis/blue_test.R ./output/blue/bluefits.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

$(blue_figs) : ./analysis/blue_plots.R ./output/blue/blue_df.csv
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

## Hypothesis plot ----

.PHONY : hyp
hyp : ./output/hyp/ternary.pdf

./output/hyp/ternary.pdf : ./analysis/hypothesis_plot.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout
