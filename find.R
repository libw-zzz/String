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
