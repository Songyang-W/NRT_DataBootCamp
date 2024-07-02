#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 13:30:33 2024

@author: songyangwang
"""

def overlay_cluster(data, clusters):
    plt.scatter(data[:, 0], data[:, 1], c=clusters, cmap='viridis')
    plt.show()

def plot_cluster(data, condition, average_plot, exclude, flag, file_name):
    # Implement the function based on the MATLAB code functionality
    pass