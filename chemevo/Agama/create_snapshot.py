import agama

mass_unit = (1.0/4.3)*(10.0**(6.0))

agama.setUnits(mass=mass_unit, length=1, velocity=1)

pot = agama.Potential(type='Spheroid', gamma=1.0, beta=3.1, scaleRadius=2.5, outerCutoffRadius=15.0)
df = agama.DistributionFunction(type='QuasiSpherical',potential=pot)
model = agama.GalaxyModel(pot,df)
M = model.sample(10000)

print(M[0][9999,0])
agama.writeSnapshot('test_snapshot.snp',M,'n')


