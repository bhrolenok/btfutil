# btf_img_to_pos.py
import btfutil
import sys

def main(btf_path,pixels_per_m,stamps_per_s):
	btf = btfutil.BTF(btf_path)
	btfutil.compute_img2pos(btf,pixels_per_m)
	btfutil.compute_ts2clock(btf,stamps_per_s)
	btf.save_to_dir(btf_path,columns=['xpos','ypos','clocktime'],overwrite=True)

if __name__ == '__main__':
	main(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]))