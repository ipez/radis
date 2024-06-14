import numpy as np

# 读取数据
input_file = 'radiance_noslit_3000.txt'
output_file = input_file[:-4]+'_interp.txt'

data = np.loadtxt(input_file, skiprows=1)  # 跳过第一行标题
wavelengths = data[:, 0]
radiances = data[:, 1]

# 创建插值后的波长
interp_wavelengths = np.arange(wavelengths[0], wavelengths[-1], 1)

# 保持已有数据点不变
interp_radiances = np.zeros_like(interp_wavelengths)
for wl, rad in zip(wavelengths, radiances):
    closest_idx = np.abs(interp_wavelengths - wl).argmin()
    interp_radiances[closest_idx] = rad

# 保存到新文件
with open(output_file, 'w') as f:
    f.write("# Wavelength [air] (nm)\tradiance_noslit (cm-1)\n")
    for wl, rad in zip(interp_wavelengths, interp_radiances):
        f.write(f"{wl:.9e}\t{rad:.9e}\n")

print(f"插值完成，结果已保存到 {output_file}")

from radis import Spectrum
s = Spectrum.from_array(interp_wavelengths, interp_radiances, 'radiance_noslit', wunit='nm', Iunit="cm-1")
s.apply_slit(5, 'nm', shape="gaussian")
filename = input_file.replace('noslit','slit')
s.savetxt(filename, 'radiance', wunit='nm', Iunit='cm-1')
print(f"狭缝卷积，结果已保存到 {filename}")

s.plot('radiance')
