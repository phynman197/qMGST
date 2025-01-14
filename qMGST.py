#!/usr/bin/env python
# Migrated to Ocean SDK from old SAPI code

import numpy as np
import sys
import traceback
import dimod
from dwave.system import DWaveSampler, FixedEmbeddingComposite, EmbeddingComposite

from dwave.cloud import Client
import time
# 使用你的 API token
client = Client.from_config(token='')

# 列出所有可用的求解器
solvers = client.get_solvers()

print("Available solvers:")
for solver in solvers:
    print(f"- {solver.id}")


# 原有代码参数
s, s2 = 1.0, 1.0
print('Embed scale=', s, s2)
assert 0 <= s <= 1

# 读取输入
# 第一行：n
# line = sys.stdin.readline().strip().split()
# n = int(line[0])
# n=45

# # 读取 QUBO 矩阵的上三角部分 (i<j)
# Q = {}
# for i in range(n):
#     line = sys.stdin.readline().strip().split()
#     for j in range(n):
#         t = float(line[j])
#         if j > i and t != 0:
#             Q[(i, j)] = t

# 定义文件路径
# 4，6，9，12
file_path = '/Users/xiaofanli/PycharmProjects/Dwave/12/qubo.txt'
matrix = np.loadtxt(file_path)
rows, cols = matrix.shape
n = rows

# 读取 QUBO 矩阵的上三角部分 (i < j)
Q = {}
with open(file_path, 'r') as file:
    lines = file.readlines()  # 读取所有行
    n = len(lines)  # 假设矩阵是 n x n 的，并且文件中有 n 行
    for i in range(n):
        line = lines[i].strip().split()  # 处理每一行
        for j in range(n):
            t = float(line[j])
            if j > i and t != 0:
                Q[(i, j)] = t

# QUBO 转为 Ising 并计算最大值用于归一化
# 注意：Ocean中dimod已内置转换工具，但这里按原来逻辑手动操作
(H, J, ising_offset) = dimod.utilities.qubo_to_ising(Q)

maxH = 0.0
if len(H):
    maxH = max(abs(min(H)), abs(max(H)))
maxJ = 0.0
if len(J.values()) > 0:
    maxJ = max(abs(min(J.values())), abs(max(J.values())))

maxV = max(maxH, maxJ) if maxH or maxJ else 1.0
if maxV == 0:
    maxV = 1.0

# 构造 BQM
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# 缩放 BQM 与原代码一致
# bqm.scale(s2 / maxV, inplace=True)
bqm.scale(s2 / maxV)

# # 读取嵌入
# embedding = eval(sys.stdin.readline())
# print('embedding=', embedding)
# qubits = sum(len(embed) for embed in embedding.values())
# print('Physical qubits used=%s' % qubits)

# 使用提供的token连接D-Wave云
# 请根据需要修改endpoint或solver名称
# 若不在环境中配置endpoint和solver，可在此显式指定，如：
# DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi', token='您的token', solver={'name': 'DW_2000Q_6'})
try:
    dwave_sampler = DWaveSampler(
        token='',
        # 若需指定endpoint和solver，请取消注释并根据实际情况填写：
        endpoint='https://cloud.dwavesys.com/sapi',
        # solver={'name': 'hybrid_binary_quadratic_model_version2'}
        solver={'name': 'Advantage2_prototype2.6'}
    )
except:
    print('Error creating DWaveSampler')
    traceback.print_exc()
    sys.exit(1)

# 使用固定嵌入
# sampler = FixedEmbeddingComposite(dwave_sampler, embedding)
sampler = EmbeddingComposite(dwave_sampler)

# 原代码参数设置
annealT, progT, readT = 160, 100, 100
print('annealT=', annealT, 'progT=', progT, 'readT=', readT)
start_time = time.time()
# 在 Ocean 中对BQM直接求解，无需手动 solve_ising
# 原代码中使用了postprocess='optimization'，Ocean中不支持该参数，
# 可考虑手动后处理，但这里先不使用后优化。
result = sampler.sample(
    bqm,
    num_reads=1000,
    annealing_time=annealT,
    programming_thermalization=progT,
    readout_thermalization=readT,
    # 对于链强度、spin_reversal_transform等参数可根据需要添加，如：
    # chain_strength=1.0, # 若需链强度，根据原代码embedding逻辑选择适当值
    # num_spin_reversal_transforms=20  # 原代码使用了20个SRT
    # num_spin_reversal_transforms=20
)

# print('result:', result)

# 使用 FixedEmbeddingComposite 后结果已是解嵌过的，不需要 unembed_answer
# 直接读取样本
samples = result.samples()
newresult = [list(sample.values()) for sample in samples]
# print('newresult:', newresult)
end_time = time.time()
running_time = end_time - start_time
print(f"程序运行时间: {running_time:.4f} 秒")

# 假设 result 是从 sampler.sample() 方法返回的结果对象
print(result.info)
# print("Total real time:", result.info['timing']['total_real_time'], "microseconds")
# print("QPU access time:", result.info['timing']['qpu_access_time'], "microseconds")
