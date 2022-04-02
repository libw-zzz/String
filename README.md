# String
Find character
## 匹配序列
1. 拆分原氨基酸序列为N段子序列 
假定原序列共有L个氨基酸，分别以子序列的开始位置(m,m最小取值为1，最大取值为原始氨基酸序列长度(L)减去最小匹配模式的长度(l,本例中为21))，和结束位置(n,n的最小取值为所匹配氨基酸的长度(l )，最大值为原始氨基酸的长度(L))对原序列进行拆分(m < n).
2. 枚举所有可能的匹配模式
在本例中，所要匹配的目标序列为XXXXKXXXN{}XXXXXN{}XKXXXR，其中两个{}内的氨基酸为任意0个或多个氨基酸，因此需要对其组合模式进行枚举，但是二者所包含的氨基酸总数不得不得超过原始氨基酸序列和所匹配序列的最短长度。
3. 依据匹配模式来挑选原氨基酸序列所拆分的子序列
