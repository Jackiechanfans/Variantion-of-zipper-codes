# A class of variant of zipper codes
A C++ implementation for simulating the encoding and decoding of the variant of zipper codes, such as multi-chain tiled diagonal zipper codes (**MC-TDZCs**) and generalized integrated interleaved zipper (**GII-zipper**) codes, and simulating the anchor decoding (**AD**) for the zipper codes.

## Reference
More details about **zipper codes**, **MC-TDZCs**, **GII-zipper codes**, and **AD** can be found in the following papers.
- Y. Sukmadji, U. Martínez-Peñas, and F. R. Kschischang, “[Zipper codes](https://ieeexplore.ieee.org/document/9843869),” J. Lightw. Technol., vol. 40, no. 19, pp. 6397-6407, Jul. 2022.
- X. Zhao, S. Zhao, and X. Ma, “[A class of tiled diagonal zipper codes with multiple chains](https://ieeexplore.ieee.org/document/9802520),” IEEE Trans. Commun., vol. 70, no. 8, pp. 5004-5017, Aug. 2022.
- X. Zhao, S. Zhao, and Z. Li, “[Enhanced anchor decoder for staircase codes with hard reliability scores](https://ieeexplore.ieee.org/document/9887957),” IEEE Commun. Lett., pp. 1–1, 2022.
- C. Häger and H. Pfister, “[Approaching miscorrection-free performance of product codes with anchor decoding](https://ieeexplore.ieee.org/document/8316914),” IEEE Tran. Commun., vol. 66, no. 7, pp. 2797-2808, Jul. 2018.

## Features
The program simulates the variant of **tiled diagonal zipper codes** (TDZCs) with **iterative bounded distance decoding** (iBDD) or **anchor decoding** (AD) over the **binary symmetric channel**,  where the size of 'tile' can be configured by the user.

The program can achieve the following simulations:
- TDZC with iBDD or AD
- MC-TDZC with iBDD or AD
- GII(4,3)-TDZC with AD

The program will run correctly under the following assumptions:
- Only BCH codes are used for all simulations
- The same interleaver map is used for each TDZC
- The same conponent code is used for each row of a TDZC
- The same buffer width for each TDZC, i.e., the same code length
- For the simulation of GII-zipper codes, only the case of the GII(4,3)-TDZC code does work

Some flexibilities are allowed:
- The code length and error correcting capacity of BCH code are allowed to be designed
- The decoding window size, the output size, and the 'tile' size are allowed to be configured by the user

If you want to get the program of other zipper schemes, like braided codes, or design the interleaver map of zipper codes, please refer to the [program](https://www.comm.utoronto.ca/frank/ZipperSim/index.html) provided by the authors of [zipper codes](https://ieeexplore.ieee.org/document/9843869).

## Configuration
There are two configuration files: 
- `GIISetting.txt`
- `Zippersetting.txt`

Make sure these two files are in the **same directory** as the source files.

### 1) TDZC with iBDD or AD
If you want to simulate **TDZC with iBDD**, please configure **GIISetting.txt** and **Zippersetting.txt**.

In `Zippersetting.txt`, make sure that:
1. The flag of "Runing_anchor_decoding_option" must be set to 0.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will not affect the simulation results.
3. The value of "Number_of_zipper_codes" must be same with the value of "Number_of_zippper_codes" in GIISetting.txt.
4. The flag of "Multichain_couple" must be set to 0.

In `GIISetting.txt`, make sure that:
1. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_zippper_codes" in Zippersetting.txt.
2. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_nested_interleaves".
3. $n < 2^q$.
4. $n$ is divisible by $w$.
5. The number of value for "The_error_correcting_capability" must be $v+1$ and the first value of "The_error_correcting_capability" must be 0. The remaining values represent the error correcting capabilities for different zipper codes.

If you want to simulate **TDZC with AD**, please configure **GIISetting.txt** and **Zippersetting.txt**.

In `Zippersetting.txt`, make sure that:
1. The flag of "Runing_anchor_decoding_option" must be set to 1.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will affect the simulation results.

Other configurations are similar to that of the cases of TDZC with iBDD.

### 2) MC-TDZC with iBDD or AD
If you want to simulate **MC-TDZC with iBDD**, please configure **GIISetting.txt** and **Zippersetting.txt**.

In `Zippersetting.txt`, make sure that:
1. The flag of "Runing_anchor_decoding_option" must be set to 0.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will not affect the simulation result.
3. The value of "Number_of_zipper_codes" must be same with the value of "Number_of_zippper_codes" in GIISetting.txt.
4. The flag of "Multichain_couple" must be set to 1.

In `GIISetting.txt`, make sure that:
1. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_zippper_codes" in Zippersetting.txt.
2. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_nested_interleaves".
3. $n < 2^q$.
4. $n$ is divisible by $w$.
5. The number of value for "The_error_correcting_capability" must be v+1 and the first value of "The_error_correcting_capability" must be 0. The remaining values represent the error correcting capabilities for different zipper codes.

If you want to simulate **MC-TDZC with AD**, please configure **GIISetting.txt** and **Zippersetting.txt**.

In `Zippersetting.txt`, make sure that:
1. The flag of "Runing_anchor_decoding_option" must be set to 1.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will affect the simulation results.

Other configurations are similar to that of the cases of MC-TDZC with iBDD.

### 3) GII(4,3)-TDZC with AD
If you want to simulate **GII(4,3)-TDZC with AD**, please configure **GIISetting.txt** and **Zippersetting.txt**.

In `Zippersetting.txt`, make sure that:
1. The flag of "Runing_anchor_decoding_option" must be set to 1.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will affect the simulation result.

In `GIISetting.txt`, make sure that:
1. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_zippper_codes" in Zippersetting.txt.
2. $n < 2^q$.
3. $n$ is divisible by $w$.

## License
The codes in this repository are under the MIT license as specified by the LICENSE file.

## Notes
If you find these codes helpful, please cite these papers:
- Y. Sukmadji, U. Martínez-Peñas, and F. R. Kschischang, “[Zipper codes](https://ieeexplore.ieee.org/document/9843869),” J. Lightw. Technol., vol. 40, no. 19, pp. 6397-6407, Jul. 2022.
- X. Zhao, S. Zhao, and X. Ma, “[A class of tiled diagonal zipper codes with multiple chains](https://ieeexplore.ieee.org/document/9802520),” IEEE Trans. Commun., vol. 70, no. 8, pp. 5004-5017, Aug. 2022.
- X. Zhao, S. Zhao, and Z. Li, “[Enhanced anchor decoder for staircase codes with hard reliability scores](https://ieeexplore.ieee.org/document/9887957),” IEEE Commun. Lett., pp. 1–1, 2022.
- C. Häger and H. Pfister, “[Approaching miscorrection-free performance of product codes with anchor decoding](https://ieeexplore.ieee.org/document/8316914),” IEEE Tran. Commun., vol. 66, no. 7, pp. 2797-2808, Jul. 2018.

If you have any questions and suggestions, please contact me at xz_xinwei@163.com.
