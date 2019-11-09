import clamato_read
import cuts
import milkyway
import prep
import bootstrap

#the high_z bin!!

spec, red, coord = clamato_read.xspec([2.5,2.75]) 

spec, red, coord = cuts.lyman_limit(spec, red, coord, .2) 

spec, red, coord = cuts.s2n(spec, red, coord, 1.5) 

spec = milkyway.extinction(spec, red, coord)

rest_spec = prep.clean(spec, red)

print(rest_spec.nspec)

bootstrap.stacks(rest_spec,red,50,"hiz")

