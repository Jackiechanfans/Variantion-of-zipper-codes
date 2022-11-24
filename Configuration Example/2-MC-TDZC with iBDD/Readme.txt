In this case, please ensure that

Zippersetting.txt
1. The flag of "Runing_anchor_decoding_option" must be set to 0.
2. The number of "Anchor_decoding_threshold_value" must be consistent with the value of "Number_of_zipper_codes", where these values of "Anchor_decoding_threshold_value" will not affect the simulation result.
3. The value of "Number_of_zipper_codes" must be same with the value of "Number_of_zippper_codes" in GIISetting.txt.
4. The flag of "Multichain_couple" must be set to 1.

GIISetting.txt
1. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_zippper_codes" in Zippersetting.txt.
2. The value of "Number_of_zippper_codes" must be same with the value of "Number_of_nested_interleaves".
3. n < 2^q.
4. n is divisible by w.
5. The number of value for "The_error_correcting_capability" must be v+1 and the first value of "The_error_correcting_capability" must be zero. The remaining values represent the error correcting capabilities for different zipper codes.