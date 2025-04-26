#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Improved QKD Network Implementation with Wavelength Division Multiplexing
"""

import SingleLink as sl
import numpy as np
from typing import Tuple, List, Dict

class QKDNetwork:
    def __init__(self, num_users: int = 6):
        self.num_users = num_users
        self.fiber_length = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]  # 6 users connected to QNSP
        
        # 预定义的波长分配矩阵
        self.wavelength_matrix1 = np.array([
            [4, -4, -2, 4, -4, -2],
            [3, 2, -3, 3, 2, -3],
            [1, 1, 1, -1, -1, -1]
        ])
        
        self.wavelength_matrix2 = np.array([
            [5, -5, -3, 5, -5, -3],
            [4, 3, -4, 4, 3, -4],
            [2, 1, 1, -2, -1, -1]
        ])
        
        # 当前使用的波长分配
        self.current_wavelength_allocation = self.wavelength_matrix1
        
        # 网络性能参数
        self.performance_metrics = {
            'average_key_rate': 0,
            'min_key_rate': 0,
            'channel_efficiency': 0,
            'noise_level': 0
        }

    def ideal_single_counts(self, wavelength_allocation: np.ndarray, user: int,
                          fiber_length: List[float]) -> Tuple[int, float]:
        """
        计算理想情况下的单光子计数
        """
        ideal_single_link = np.zeros(self.num_users, dtype=int)
        pgr_per_channel = 1E6 / np.max(wavelength_allocation)
        
        singlecounts = sl.Calc_QKD_Parameters_For_Single_Link(
            PairGenerationRate_per_second=pgr_per_channel,
            FibreLength_in_km=[fiber_length[user], 0],
            SingleCountsFromOtherUsers_per_second=[0, 0, 0, 0]
        )
        
        ideal_single_link = singlecounts['Singles rate for Alice'].HA / 0.7
        return ideal_single_link, pgr_per_channel

    def wavelength_demultiplexing(self, node_index: int, 
                                wavelength_channels: List[int]) -> Dict:
        """
        实现接收端的波长去复用
        """
        demuxed_channels = {}
        for wavelength in wavelength_channels:
            # 为每个波长计算最优的符合窗口
            window = self._calculate_optimal_window(wavelength)
            # 为每个波长分配探测器
            detector = self._assign_detector(wavelength)
            
            demuxed_channels[wavelength] = {
                'window': window,
                'detector': detector
            }
        return demuxed_channels

    def _calculate_optimal_window(self, wavelength: int) -> float:
        """
        计算给定波长的最优符合窗口
        """
        base_window = 150  # ps
        # 根据波长调整窗口大小
        dispersion_factor = abs(wavelength) * 0.1
        return base_window + dispersion_factor

    def _assign_detector(self, wavelength: int) -> int:
        """
        为波长分配探测器
        """
        # 简单的轮询分配策略
        return abs(wavelength) % 2

    def num_channel_from_other_users(self, wavelength_allocation: np.ndarray,
                                   row: int, col: int) -> Tuple[List, List, int, int]:
        """
        计算来自其他用户的通道数
        """
        channels_at_nodes = []
        effective_channel = []
        
        # 收集节点上的所有通道
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if y == row or y == col:
                    channels_at_nodes.append(wavelength_allocation[x][y])
        
        # 计算有效通道
        num_eff_chan = 0
        if row > col:
            for num in channels_at_nodes[0::2]:
                if -num in channels_at_nodes:
                    effective_channel.append(num)    
                    num_eff_chan += 1
        else:
            for num in channels_at_nodes[1::2]:
                if -num in channels_at_nodes:
                    effective_channel.append(num)    
                    num_eff_chan += 1
        
        # 计算噪声通道
        num_noise_chan = 0
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if wavelength_allocation[x][y] * (-1) in effective_channel:
                    num_noise_chan += 1
        
        num_noise_chan = num_noise_chan - len(effective_channel)
        return channels_at_nodes, effective_channel, num_eff_chan, num_noise_chan

    def node_noise_counts(self, wavelength_allocation: np.ndarray) -> np.ndarray:
        """
        计算节点噪声计数
        """
        node_others = np.zeros((self.num_users, self.num_users), dtype=int)
        
        for row in range(self.num_users):
            for col in range(self.num_users):
                single_counts, _ = self.ideal_single_counts(
                    wavelength_allocation, row, self.fiber_length
                )
                _, _, _, num_noise_chan = self.num_channel_from_other_users(
                    wavelength_allocation, row, col
                )
                node_others[row][col] = single_counts * num_noise_chan
                
        return node_others

    def multi_secure_key_rate(self, wavelength_allocation: np.ndarray,
                            coincidence_window: int) -> np.ndarray:
        """
        计算多用户情况下的安全密钥生成率
        """
        secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        node_others = self.node_noise_counts(wavelength_allocation)
        
        for row in range(self.num_users):
            for col in range(self.num_users):
                if row == col:
                    continue
                    
                _, pgr_per_channel = self.ideal_single_counts(
                    wavelength_allocation, row, self.fiber_length
                )
                _, _, num_eff_chan, _ = self.num_channel_from_other_users(
                    wavelength_allocation, row, col
                )
                
                new_single_link = sl.Calc_QKD_Parameters_For_Single_Link(
                    PairGenerationRate_per_second=pgr_per_channel * num_eff_chan,
                    CoincidenceWindow_in_ps=coincidence_window,
                    FibreLength_in_km=[self.fiber_length[row], self.fiber_length[col]],
                    SingleCountsFromOtherUsers_per_second=[
                        node_others[row][col], node_others[row][col],
                        node_others[col][row], node_others[col][row]
                    ]
                )
                
                secure_key_rate_matrix[row][col] = new_single_link['Total secure key rate']
                
        return secure_key_rate_matrix

    def optimize_coincidence_window(self, wavelength_allocation: np.ndarray) -> Tuple[np.ndarray, int, int]:
        """
        优化符合窗口以获得最大总安全密钥率
        """
        optimized_secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        optimized_cd = [0, 0]  # [最佳总密钥率, 最佳符合窗口]
        
        for cd in range(10, 400, 1):
            skr_matrix = self.multi_secure_key_rate(wavelength_allocation, cd)
            sum_of_secure_key_rate = np.sum(skr_matrix)
            
            if sum_of_secure_key_rate > optimized_cd[0]:
                optimized_cd[0] = sum_of_secure_key_rate
                optimized_cd[1] = cd
                optimized_secure_key_rate_matrix = skr_matrix
                
        return optimized_secure_key_rate_matrix, optimized_cd[1], optimized_cd[0]

    def evaluate_network_performance(self, secure_key_rate_matrix: np.ndarray) -> Dict:
        """
        评估网络性能
        """
        non_zero_rates = secure_key_rate_matrix[secure_key_rate_matrix > 0]
        
        self.performance_metrics = {
            'average_key_rate': np.mean(non_zero_rates),
            'min_key_rate': np.min(non_zero_rates),
            'channel_efficiency': len(non_zero_rates) / (self.num_users * (self.num_users - 1)),
            'noise_level': self._calculate_noise_level()
        }
        
        return self.performance_metrics

    def _calculate_noise_level(self) -> float:
        """
        计算网络噪声水平
        """
        node_others = self.node_noise_counts(self.current_wavelength_allocation)
        return np.mean(node_others)

    def dynamic_wavelength_allocation(self) -> np.ndarray:
        """
        根据网络状态动态调整波长分配
        """
        current_performance = self.evaluate_network_performance(
            self.multi_secure_key_rate(self.current_wavelength_allocation, 150)
        )
        
        # 测试不同的波长分配方案
        test_allocations = [self.wavelength_matrix1, self.wavelength_matrix2]
        best_allocation = self.current_wavelength_allocation
        best_performance = current_performance['average_key_rate']
        
        for allocation in test_allocations:
            test_skr = self.multi_secure_key_rate(allocation, 150)
            test_performance = np.mean(test_skr[test_skr > 0])
            
            if test_performance > best_performance:
                best_performance = test_performance
                best_allocation = allocation
        
        return best_allocation

def main():
    # 创建网络实例
    network = QKDNetwork(num_users=6)
    
    # 优化当前波长分配的性能
    opt_skr_matrix, opt_cd, opt_tskr = network.optimize_coincidence_window(
        network.current_wavelength_allocation
    )
    
    # 评估网络性能
    performance = network.evaluate_network_performance(opt_skr_matrix)
    
    # 打印结果
    print('Optimized coincidence window:', opt_cd)
    print('Optimized secure key rate matrix:\n', opt_skr_matrix)
    print('Total secure key rate:', opt_tskr)
    print('Network performance metrics:', performance)
    
    # 尝试动态波长分配
    new_allocation = network.dynamic_wavelength_allocation()
    if not np.array_equal(new_allocation, network.current_wavelength_allocation):
        print('New wavelength allocation suggested:\n', new_allocation)

if __name__ == "__main__":
    main()