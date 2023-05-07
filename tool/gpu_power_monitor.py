"""
================================================================================
@file:        gpu_power_monitor.py
@brief:       Monitors and prints real-time GPU power usage with timestamps.
@author:      Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
              Oh-Kyoung Kwon (okkwon@kisti.re.kr), Korea Institute of Science and Technology Information
              Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
@date:        May 2023
@version:     2.0
@copyright:   Copyright (c) 2023 Mingyu Yang, Oh-Kyoung Kwon, Jung-Il Choi. All rights reserved.
@license:     This project is released under the terms of the MIT License (see LICENSE file).
================================================================================
"""

import pynvml
import matplotlib.pyplot as plt
import time
import numpy as np
import asyncio
from datetime import datetime, timedelta

# Initialize NVIDIA Management Library (NVML)
pynvml.nvmlInit()

# Get handle for the first available GPU
handle = pynvml.nvmlDeviceGetHandleByIndex(0)

# Set measurement interval (in seconds)
measurement_interval = 0.01
mW_to_W = 1e3  # Conversion factor for milliwatts to watts

# Continuously monitor and print GPU power usage
while True:
    now = datetime.now()
    current_datetime = now.strftime("%Y%m%d %H:%M:%S.%f")
    sec = 0.01
    measure = pynvml.nvmlDeviceGetPowerUsage(handle) / mW_to_W  # Get power usage and convert to watts
    time.sleep(sec)
    print("time : ", current_datetime, measure)  # Print timestamp and power usage
