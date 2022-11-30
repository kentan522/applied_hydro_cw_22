from WaveCase import WaveCase

# Defining wave cases
wave_cases = [
    (27.5, 14),
    (24.7, 12),
    (19.8, 10),
]

# Instantiating wave cases
stored_cases = {}
for idx, (h_max, T) in enumerate(wave_cases):
    print('Case ', (idx + 1), ':')
    stored_cases[idx] = WaveCase(h_max, T)
    print('')


print('end of program')




