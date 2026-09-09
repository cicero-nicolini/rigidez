import rigidez

conc1 = rigidez.Concreto(fck=20)

conc1.fck = 90

print(conc1.ηc)
print(conc1.fck)