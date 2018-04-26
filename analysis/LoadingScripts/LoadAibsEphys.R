## load aibs ephys data


aibsEphys <- read.csv("data-raw/aibs_ephys/aibs_aggregated_ephys_v6.csv")
aibsEphys$transgene = aibsEphys$transgene %>% make.names
aibsEphys$transgene  = factor(aibsEphys$transgene)

aibsEphys$NeuronType = aibsEphys$transgene
aibsEphys = rename(aibsEphys, sample.name = transgene)
aibsEphys$rmp = aibsEphys$vrest
aibsEphys$rin = aibsEphys$input_resistance_mohm
aibsEphys$tau = aibsEphys$tau
aibsEphys$apthr = aibsEphys$threshold_v_long_square
aibsEphys$apamp = aibsEphys$peak_v_long_square - aibsEphys$threshold_v_long_square
aibsEphys$ahpamp = aibsEphys$threshold_v_long_square - aibsEphys$fast_trough_v_long_square
aibsEphys$aphw = aibsEphys$rheo_first_spike_hw * 1000
#aibsEphys$aphw = aibsEphys$upstroke_downstroke_ratio_long_square
aibsEphys$rheo = aibsEphys$threshold_i_long_square
aibsEphys$cap = (aibsEphys$tau / aibsEphys$rin) * 1000
aibsEphys$maxfreq = aibsEphys$max_rate_long_square#1 / (aibsEphys$avg_isi / 1000)
aibsEphys$adratio = aibsEphys$max_rate_first_mean_adratio# 1 / ((aibsEphys$adaptation * 10) +1)
#1/((1 + (aibsEphys$adaptation * 10) # note that aibs is actually calculating the adaptation index, not ratio

aibsEphys = merge(aibsEphys, aibsColors, by = "sample.name")

saveRDS(aibsEphys, file = 'data/aibs_ephys.rda')
