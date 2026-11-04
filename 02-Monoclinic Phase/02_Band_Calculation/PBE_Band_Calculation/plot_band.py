import matplotlib.pyplot as plt
import numpy as np

# 1. KLABELS dosyasından yüksek simetri noktalarını oku
k_ticks = []
k_labels = []
with open("KLABELS", "r") as f:
  lines = f.readlines()[1:]  # Başlık satırını atla
  for line in lines:
    parts = line.strip().split()
    if len(parts) >= 2 and not line.startswith("*"):
      try:
        pos = float(parts[-1])
        label = parts[0]
        if label == "GAMMA" or label == "G":
          label = r"$\Gamma$"
        k_ticks.append(pos)
        k_labels.append(label)
      except ValueError:
        continue

# 2. REFORMATTED_BAND.dat dosyasını oku (1. sütun: k-mesafesi, diğer sütunlar: bantlar)
data = np.loadtxt("REFORMATTED_BAND.dat")
k_dist = data[:, 0]
bands = data[:, 1:]

# 3. Grafiği oluştur
fig, ax = plt.subplots(figsize=(7, 5), dpi=300)

for i in range(bands.shape[1]):
  ax.plot(k_dist, bands[:, i], color="#1f77b4", lw=1.2)

# Yüksek simetri dikey çizgileri ve Fermi seviyesi
for tick in k_ticks:
  ax.axvline(x=tick, color="gray", linestyle="--", lw=0.7, alpha=0.7)
ax.axhline(y=0.0, color="red", linestyle="--", lw=1.0, label="$E_F$")

ax.set_xticks(k_ticks)
ax.set_xticklabels(k_labels, fontsize=11)
ax.set_xlim(k_dist.min(), k_dist.max())
ax.set_ylim(-2.5, 2.5)  # Dar bant aralığını net görmek için [-2.5, 2.5] eV
ax.set_ylabel("Energy ($E - E_F$) [eV]", fontsize=12)
ax.set_title(r"$\beta$-Ag$_2$Te Band Structure (PBE+SOC+vdW)", fontsize=13)

plt.tight_layout()
plt.savefig("PBE_SOC_BandStructure.png", dpi=300)
plt.show()
