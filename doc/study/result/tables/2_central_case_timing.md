| case | implementation | mode | np | iter1_9_total_ms_mean | compute_ms_mean | comm_ms_mean | packing_ms_mean | throughput_Mcells_s | iter0_to_iter1_9_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| np2_128x128x4096_device | cuda-cxx | device | 2 | 4.225 | 4.114 | 0.052 | 0.060 | 15884.7 | 10.7 |
| np4_128x128x4096_device | cuda-cxx | device | 4 | 2.221 | 2.058 | 0.070 | 0.097 | 30214.4 | 53.0 |
| np8_128x128x4096_device | cuda-cxx | device | 8 | 1.340 | 1.042 | 0.147 | 0.171 | 50095.9 | 160.1 |
| np2_128x128x4096_host | cuda-cxx | host | 2 | 4.600 | 4.115 | 0.426 | 0.060 | 14590.3 | 6.9 |
| np4_128x128x4096_host | cuda-cxx | host | 4 | 2.438 | 2.058 | 0.288 | 0.096 | 27531.7 | 28.8 |
| np8_128x128x4096_host | cuda-cxx | host | 8 | 1.650 | 1.030 | 0.458 | 0.174 | 40677.0 | 89.8 |
