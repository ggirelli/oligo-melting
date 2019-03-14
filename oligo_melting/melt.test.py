from melt import *

R = 1.987 / 1000    # kcal / (K mol)
oligo = 2.5e-6

dimers = ["AT", "TG", "GT", "TC", "CA", "AA", "AT", "TG", "GT", "TC",
	"CA", "AA", "AT", "TG", "GT", "TC", "CA", "AA", "AT", "TG", "GT", "TC", "CA"]

nnet = NN_TABLES['DNA:DNA']

h0 = np.sum([nnet.dH0[x] for x in dimers])
s0 = np.sum([nnet.dS0[x] for x in dimers])
s0 /= 1e3

tm = h0 / (s0 + R * np.log(oligo)) - 273.15

print(tm)