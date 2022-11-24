In this case, please ensure that

Zippersetting.txt
1. The flag of "Runing_anchor_decoding_option" must be set to 1.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will affect the simulation result.

GIISetting.txt
1. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_zipper_codes" in Zippersetting.txt.
2. n < 2^q.

Note that this program only support the case of GII(4,3)-TDZC code. If you want to get other cases, like GII(3,1)-TDZC code, please contact xz_xinwei@163.com.