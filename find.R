### 匹配序列
### 拆分原氨基酸序列为N段子序列 假定原序列共有L个氨基酸，分别以子序列的开始位置(m,m最小取值为1，最大取值为原始氨基酸序列长度(L)减去最小匹配模式的长度(l,本例中为21))，和结束位置(n,n的最小取值为所匹配氨基酸的长度(l )，最大值为原始氨基酸的长度(L))对原序列进行拆分(m < n).
### 枚举所有可能的匹配模式 在本例中，所要匹配的目标序列为XXXXKXXXN{}XXXXXN{}XKXXXR，其中两个{}内的氨基酸为任意0个或多个氨基酸，因此需要对其组合模式进行枚举，但是二者所包含的氨基酸总数不得不得超过原始氨基酸序列和所匹配序列的最短长度。
### 依据匹配模式来挑选原氨基酸序列所拆分的子序列

### XXXXKXXXN...XXXXXN...XKXXXR
DNA_seq <- "MVLIESSESEDEILIKNEPKPASSSSPQPTETKQVDGDDSDGFETASEREISDEEGEEDGTKNDAVTSQEEPQHSEKKEEQIELMSEGEAIVDDGSNKEKALAEANEAKAEGNKLFVNGLYEEALSKYAFALELVQELPESIELRSICYLNRGVCFLKLGKCEETIKECTKALELNPTYNKALVRRAEAHEKLEHFEDAVTDLKKILELDPSNDQARKGIRRLEPLAAEKREKMKEEAITKLKEMGNSILGRFGMSVDNFKAVKDPNTGSYSLSFQN"

DNA_seq_length <- nchar(DNA_seq)


###拆分原始DNA序列为N段子序列
DNA_seq_seg <- array(NA, dim = c(DNA_seq_length, DNA_seq_length)) %>% 
  as.data.frame() %>% 
  mutate(Start = 1:DNA_seq_length)

for(i in 1:(DNA_seq_length-20)){
  for(j in 21:DNA_seq_length){
    DNA_seq_seg[i, j] <- str_sub(DNA_seq, i, j)
  }
}

DNA_seq_seg %>% 
  as.data.frame() %>% 
  pivot_longer(1:DNA_seq_length, values_to = "dna_seg", names_to = "End") %>% 
  mutate_at(vars(End), ~gsub("V", "", .)) %>% 
  mutate(Seg_length = nchar(dna_seg)) %>% 
  mutate_at(vars(Seg_length), ~ifelse(.<21, NA, .)) %>% 
  na.omit() %>% 
  unite(dna_loc_seg, Start, End, dna_seg) %>% 
  group_by(Seg_length) %>% 
  nest() -> DNA_seg_nest


### 枚举匹配模式
str_leng <- DNA_seq_length - 15
pattern_n <- array(NA, dim = c(str_leng - 4, str_leng)) %>% 
  as.data.frame() %>% mutate(X = 5:str_leng)
for(i in 5:(str_leng)){
  for(j in 1:str_leng){
    
    if(i + j <= str_leng){
    
    pattern_n[i-4, j] <- paste("\\d+_\\d+_[GAVLIPFYWSTCMNQDEKRH]{4}[K]",
                 "[GAVLIPFYWSTCMNQDEKRH]{3}[N]",
                 "[GAVLIPFYWSTCMNQDEKRH]{", i, "}[N]",
                 "[GAVLIPFYWSTCMNQDEKRH]{", j, "}[K]",
                 "[GAVLIPFYWSTCMNQDEKRH]{3}[R]", sep = "")
    } else {
    pattern_n[i-4, j] <- NA
    }
  }
}

pattern_n %>% 
  pivot_longer(cols = 1:str_leng, names_to = "V2", values_to = "Pattern") %>%
  mutate(Y = gsub("V", "", V2)) %>% 
  na.omit() %>% 
  mutate_at(vars(Y), ~as.numeric(.)) %>% 
  mutate(Sum = X + Y) %>% 
  dplyr::filter(Sum <= str_leng) -> pattern2


pattern_dna <- pattern2$Pattern

map(pattern_dna, function(var){
  res1 <- str_match_all(var, "(?<=\\{)\\d+(?=\\})") %>% unlist() %>% as.numeric(.) %>% sum() + 5
}) %>% unlist() -> res2

pattern_dna_catolog <- cbind(res2, pattern_dna) %>% as.data.frame()

pattern_dna_catolog %>% 
  group_by(res2) %>% 
  nest() -> pattern_nest

### 进行匹配
match <- function(DNA_seq, pattern){
    str_match(DNA_seq, pattern) %>% as.data.frame() %>% na.omit()
}


x <- DNA_seg_nest$data
y <- pattern_nest$data

map2(x, y, function(x, y){
  
  res <- vector("list", length(y))

  for(i in 1:length(y$pattern_dna)){
    
    res[[i]] <- match(x$dna_loc_seg, y$pattern_dna[[i]])
    
  }
  
  res %>% ldply() %>% as.data.frame() %>% na.omit() %>% unique()

}) %>% ldply() %>% 
  separate(V1, c("Start", "End", "DNA_Segment"), sep = "_") -> DNA_Seg_TAR
