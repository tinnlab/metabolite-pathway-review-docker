'''

迭代边

作用：比如将最小生成树映射到原始图

'''
def itera(a):
    m = []
    for i in range(len(a) - 1):
        m.append((a[i],a[i+1]))
    return m