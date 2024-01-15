nc = 6
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

.PHONY : all
all : sims blue

## Simulations ----

.PHONY : sims
sims : ./output/sims/gsims.csv ./output/sims/glsims.csv ./output/sims/g_altsims.csv ./output/sims/gl_altsims.csv

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
