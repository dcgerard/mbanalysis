nc = 6
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

alt_sim_figs = ./output/sims/plots/t1e.tex \
               ./output/sims/plots/alt_p_box.pdf

.PHONY : all
all : sims blue

## Simulations ----

.PHONY : sims
sims : $(null_sim_figs) $(alt_sim_figs)

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

## Blueberries ----

.PHONY : blue
blue : ./output/blue/blue_df.csv

./output/blue/bluefits.RDS : ./analysis/blue_up.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

./output/blue/blue_df.csv : ./analysis/blue_test.R ./output/blue/bluefits.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout
