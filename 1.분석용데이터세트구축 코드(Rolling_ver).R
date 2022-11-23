
require(haven)
require(tidyverse)
require(tibble)
require(stringr)
require(fastDummies)
require(comprehenr) # List comprehension in R
require(glue) # f-string in R


options(tibble.width = Inf)
options(scipen = 999)
options(dplyr.summarise.inform=F)

file_path <-  'D:\\Cohort_DB\\'
save_path <- 'D:\\Sangwoo(D)\\NHIS\\Datasets_for_analysis\\DB1.0\\1.Rolling analysis datasets\\'

base_year = 2006
track_windows = 5 # 미래 발병 추적기간
past_windows = 5 # 과거 기질환 추적기간



# 약제코드 로드
dm_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\anti_dm_all.csv') %>% 
  rename(gnl_nm_cd = `Anti-diabeticmedicaiton(all)`) %>% 
  .$gnl_nm_cd
hprts_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\anti_hprts.csv') %>% 
  rename(gnl_nm_cd = `Anti-hypertensiveDrugs`) %>% 
  .$gnl_nm_cd
asthma_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\asthma.csv') %>% 
  .$gnl_nm_cd
copd_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\copd.csv') %>% 
  .$gnl_nm_cd
ra_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\ra.csv') %>% 
  .$gnl_nm_cd
statin_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\statin.csv') %>% 
  rename(gnl_nm_cd = `Statin Drugs`) %>% 
  .$gnl_nm_cd





##### 1. 설명변수 전처리


#### 1.1. 자격데이터 전처리

## 전체 자격 데이터 로드
p_table_labels <- list.files(path=file_path, pattern = 'p_20.*.dta')
p_table_list <- lapply(p_table_labels, function(x) read_dta(paste(file_path, x, sep = '')))
p_table_list2 <- lapply(p_table_list, function(x) mutate(x, DTH_CODE1 = as.character(DTH_CODE1),
                                                         DTH_CODE2 = as.character(DTH_CODE2)))
p_df <- bind_rows(p_table_list2) %>% 
  rename_all(tolower)


rm(p_table_labels)
rm(p_table_list)
rm(p_table_list2)


## Birth year 데이터 조인
birth_year <- read_dta('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\NHIS\\분석용데이터세트구축\\DB1.0\\birthyear.dta') %>% 
  rename_all(tolower)

p_df2 <- p_df %>% 
  inner_join(birth_year, by='person_id') %>% 
  select(-byr_q)




# 기준년도 자격여부
q_baseyear <- p_df2 %>% 
  filter(stnd_y == base_year) %>% 
  .$person_id



# 추적관찰기간 동안 한번도 자격DB에 없었던 person_id
no_trace_years <- p_df2 %>% 
  dummy_cols(select_columns = 'stnd_y') %>% 
  group_by(person_id, sex, byr) %>% 
  summarise_at(vars(!!sym(paste0('stnd_y_', base_year+1)):!!sym(paste0('stnd_y_', base_year+track_windows))), max) %>% 
  ungroup() %>% 
  mutate(rowsum = rowSums(select(., -c(person_id, sex, byr)))) %>% 
  filter(rowsum == 0) %>% 
  .$person_id


# 기준년도 사망여부 및 사망년도
p_death <- p_df %>% 
  mutate(dth_year = as.double(str_sub(dth_ym, 1, 4))) %>% 
  select(person_id, dth_year) %>% 
  mutate(dth_baseyear = ifelse(is.na(dth_year) | dth_year > base_year, 0, 1)) %>% 
  group_by(person_id) %>% 
  summarise(dth_baseyear = max(dth_baseyear),
            dth_year = max(dth_year, na.rm=T)) %>% 
  mutate(dth_year = na_if(dth_year, -Inf))


# 전체 종합: person_id, sex, birth year, 기준년도 자격 여부, 추적 관찰기간동안 자격DB 포함 여부, 사망
p_final <- p_df2 %>% 
  distinct(person_id, sex, byr) %>% 
  mutate(sex = ifelse(sex==2, 0, sex)) %>% 
  
  # if qualified during the base year
  mutate(q_baseyear = ifelse(person_id %in% q_baseyear, 1, 0)) %>% 
  
  # if never contained during the trace years -> 0
  # if ever contained -> 1
  mutate(c_traceyears = ifelse(person_id %in% no_trace_years, 0, 1)) %>% 
  
  # death
  left_join(p_death, by='person_id')



write_csv(p_final, paste0(save_path, '1.total_qualified_ids_', base_year, '.csv'))


rm(p_death)
rm(p_df)
rm(p_df2)



#### 1.2. 검진데이터 전처리

## 검진데이터(관찰기간) 로드 및 1차 전처리
e_df <- tibble()
for(y in base_year:base_year){
  path <- paste0('D:/Cohort_DB/e_', y, '.dta')
  if (y < 2009){
    tmp <- read_dta(path) %>% 
      rename_all(tolower) %>% 
      
      rename(year=hchk_year) %>% 
      mutate(bmi=weight/(height/100)^2) %>% 
      mutate(prote=ifelse(olig_prote_cd==1, 0, 1)) %>% 
      rename_at(vars(starts_with('fmly')), ~str_replace(., '_patien_yn', '')) %>% 
      
      mutate(smk_stat = NA,
             smk_stat = case_when(smk_stat_type_rsps_cd == 1 ~ 0,
                                  smk_stat_type_rsps_cd == 2 ~ 1,
                                  smk_stat_type_rsps_cd == 3 & dsqty_rsps_cd %in% c(1, 2) ~ 2,
                                  smk_stat_type_rsps_cd == 3 & dsqty_rsps_cd %in% c(3, 4) ~ 3)) %>% 
      mutate(drnk_dum = NA,
             drnk_dum = case_when(drnk_habit_rsps_cd <= 2 ~ 0,
                                  drnk_habit_rsps_cd > 2 ~ 1)) %>% 
      mutate(exerci = NA,
             exerci = case_when(exerci_freq_rsps_cd == 1 ~ 0,
                                exerci_freq_rsps_cd > 1 ~ 1)) %>% 
      mutate_at(vars(starts_with('fmly')), ~case_when(. == 1 ~ 0,
                                                      . == 2 ~ 1))
    
    e_df <- bind_rows(e_df, tmp)
  }
  
  if (between(y, 2009, 2017)){
    tmp <- read_dta(path) %>% 
      rename_all(tolower) %>% 
      
      rename(year=hchk_year) %>% 
      mutate(bmi=weight/(height/100)^2) %>% 
      mutate(prote=ifelse(olig_prote_cd==1, 0, 1)) %>% 
      rename_at(vars(starts_with('fmly')), ~str_replace(., '_patien_yn', '')) %>% 
      
      mutate(smk_stat = NA,
             smk_stat = case_when(smk_stat_type_rsps_cd == 1 ~ 0,
                                  smk_stat_type_rsps_cd == 2 ~ 1,
                                  smk_stat_type_rsps_cd == 3 & cur_dsqty_rsps_cd <=19 ~ 2,
                                  smk_stat_type_rsps_cd == 3 & cur_dsqty_rsps_cd >19 ~ 3)) %>% 
      mutate(drnk_dum = NA,
             drnk_dum = case_when(drnk_habit_rsps_cd==1 ~ 0,
                                  drnk_habit_rsps_cd>1 ~ 1)) %>% 
      mutate(exerci = NA,
             exerci = case_when(mov20_wek_freq_id == 1 & mov30_wek_freq_id ==1 ~ 0,
                                mov20_wek_freq_id > 1  | mov30_wek_freq_id >1  ~ 1))
    
    e_df <- bind_rows(e_df, tmp)
  }
  
  rm(tmp)
}



## 이상치 나타내는 변수 생성
e_df2 <- e_df
outlier_cols <- c('bmi', 'blds', 'tot_chole', 'sgot_ast', 'sgpt_alt', 'gamma_gtp')
mean_val <- list()
std_val <- list()

for (c in outlier_cols){
  mean_val[c] = mean(e_df[[c]], na.rm=T)
  std_val[c] = sd(e_df[[c]], na.rm=T)
}


for (c in outlier_cols){
  col_label <- paste0(c, '_outlier')
  e_df2[[col_label]] = ifelse(e_df2[[c]] <= mean_val[[c]] + 3*std_val[[c]], 0, 1)
}






## 전처리2
# 분석대상자격자 테이블 조인
p_target <- p_final %>% 
  filter(q_baseyear == 1,
         c_traceyears == 1,
         dth_baseyear == 0) %>% 
  select(person_id, sex, byr)


e_df3 <- p_target %>% 
  inner_join(e_df2, by='person_id') %>% 
  mutate(age = base_year-byr) %>% 
  select(-byr)


# 두개 년도 모두 검진을 받은 경우 최신 검진만 남기고 드랍
e_df3 <- e_df3 %>% 
  arrange(person_id, desc(year)) %>% 
  distinct(person_id, .keep_all = T)


e_df3 <- e_df3 %>% 
  # Smoke
  dummy_cols(select_columns = 'smk_stat') %>% 
  select(-c(smk_stat_0, smk_stat)) %>% 
  rename_at(vars(starts_with('smk_stat_')), ~ str_replace(., 'smk_stat_', 'smk')) %>% 
  mutate_at(vars(starts_with('smk')), ~ ifelse(is.na(.), 0, .)) %>% 
  
  # BMI boundary
  mutate(bmi_low = ifelse(bmi<18.5, bmi-18.5, 0),
         bmi_high = ifelse(bmi>23, bmi-23, 0)) %>% 
  
  # blood pressure
  mutate(bp_high_high = ifelse(bp_high >=120, bp_high -120, 0),
         bp_lwst_high = ifelse(bp_lwst >=80, bp_lwst -80, 0)) %>% 
  
  # BLDs
  mutate(blds_high = ifelse(blds>=100, blds-100, 0)) %>% 
  
  # Total Cholesterol
  mutate(tot_chole_high = ifelse(tot_chole>=200, tot_chole-200, 0)) %>% 
  
  # Hemoglobin
  mutate(hmg_low = ifelse(hmg<13 & sex==1, hmg-13, 0),
         hmg_low = ifelse(hmg<12 & sex==0, hmg-12, hmg_low)) %>% 
  
  # SGOT(AST), SGPT(ALT)
  mutate(sgot_ast_high = ifelse(sgot_ast>40, sgot_ast-40, 0),
         sgpt_alt_high = ifelse(sgpt_alt>35, sgpt_alt-35, 0)) %>% 
  
  # gamma gtp
  mutate(gamma_gtp_high = ifelse(gamma_gtp>=64 & sex==1, gamma_gtp-64, 0),
         gamma_gtp_high = ifelse(gamma_gtp>35 & sex==0, gamma_gtp-35, gamma_gtp_high)) %>% 
  
  # drnk & age
  mutate(drnk_old = ifelse(drnk_dum==1 & age>=65, 1, 0),
         drnk_young = ifelse(drnk_dum==1 & age<65, 1, 0),
         
         drnk_old = replace(drnk_old, is.na(drnk_dum), NA),
         drnk_young = replace(drnk_young, is.na(drnk_dum), NA))



# 관찰년도 중 건강검진을 받은사람의 ID
eid <- e_df3 %>% 
  .$person_id




write_csv(e_df3, paste0(save_path, '2.baseyear_checkup_', base_year, '.csv'))

rm(e_df)
rm(e_df2)



#### 1.3. 기병력 여부

## 기병력 여부
# 과거 기질환 추적기간 20T 로드
m120_df <- tibble()
for (y in max(c(2002, base_year-past_windows+1)):base_year) {
  path <- paste0(file_path, 'm120_', y, '.dta')
  
  m120 <- read_dta(path) %>% 
    rename_all(tolower)
  
  tmp <- p_final %>% 
    filter(person_id %in% eid) %>% 
    select(person_id) %>% 
    
    left_join(m120, by='person_id') %>% 
    mutate(year = as.double(str_sub(recu_fr_dt, 1, 4))) %>% 
    mutate(sub_sick = na_if(sub_sick, '')) %>% 
    mutate(main1 = str_sub(main_sick, 1, 1),
           main2 = str_sub(main_sick, 2, -1),
           sub1 = str_sub(sub_sick, 1, 1),
           sub2 = str_sub(sub_sick, 2, -1)) %>% 
    select(person_id, key_seq, year, main1, main2, sub1, sub2)
  
  
  
  m120_df <- bind_rows(m120_df, tmp)
  
  rm(tmp)
}


# 기병력 여부 변수 생성
m120_df2 <- m120_df %>% 
  # 뇌졸중
  mutate(stroke_p = ifelse((main1=='I' & str_sub(main2, 1, 2) %in% 60:69) |
                             (sub1=='I' & str_sub(sub2, 1, 2) %in% 60:69), 1, 0)) %>% 
  
  # 심근경색
  mutate(mi_p = ifelse((main1=='I' & str_sub(main2, 1, 2) %in% 21:24) |
                         (sub1=='I' & str_sub(sub2, 1, 2) %in% 21:24), 1, 0)) %>% 
  
  # 골절
  mutate(fracture_p = ifelse((main1=='S' & str_sub(main2, 1, 2) %in% c(22, 32, 72)) | (main1=='M' & str_sub(main2, 1, 2) == 84) |
                               (sub1=='S' & str_sub(sub2, 1, 2) %in% c(22, 32, 72)) | (sub1=='M' & str_sub(sub2, 1, 2) == 84), 1, 0)) %>% 
  
  # 치매
  mutate(dementia_p = ifelse((main1=='F' & str_sub(main2, 1, 2) %in% c('00', '01', '02', '03')) | (main1=='G' & str_sub(main2, 1, 2) %in% 30:32) |
                               (sub1=='F' & str_sub(sub2, 1, 2) %in% c('00', '01', '02', '03')) | (sub1=='G' & str_sub(sub2, 1, 2) %in% 30:32), 1, 0)) %>% 
  
  # 우울증
  mutate(depression_p = ifelse((main1=='F' & str_sub(main2, 1, 2) %in% 32:33) |
                                 (sub1=='F' & str_sub(sub2, 1, 2) %in% 32:33), 1, 0)) %>% 
  
  # 심부전
  mutate(chf_p = ifelse((main1=='I' & str_sub(main2, 1, 2) == 50) |
                          (sub1=='I' & str_sub(sub2, 1, 2) == 50), 1, 0)) %>% 
  
  # 폐렴
  mutate(pneumonia_p = ifelse((main1=='J' & str_sub(main2, 1, 2) %in% 12:17) |
                                (sub1=='J' & str_sub(sub2, 1, 2) %in% 12:17), 1, 0)) %>% 
  
  # 만성신장병
  mutate(ckd_p = ifelse((main1=='N' & str_sub(main2, 1, 2) == 18) |
                          (sub1=='N' & str_sub(sub2, 1, 2) == 18), 1, 0)) %>% 
  
  # 위암
  mutate(gc_p = ifelse((main1=='C' & str_sub(main2, 1, 2) == 16) |
                         (sub1=='C' & str_sub(sub2, 1, 2) == 16), 1, 0)) %>%   
  
  # 대장암
  mutate(cc_p = ifelse((main1=='C' & str_sub(main2, 1, 2) %in% 18:21) |
                         (sub1=='C' & str_sub(sub2, 1, 2) %in% 18:21), 1, 0)) %>% 
  
  # 간암
  mutate(hc_p = ifelse((main1=='C' & str_sub(main2, 1, 2) == 22) |
                         (sub1=='C' & str_sub(sub2, 1, 2) == 22), 1, 0)) %>% 
  
  # 폐암
  mutate(lc_p = ifelse((main1=='C' & str_sub(main2, 1, 2) == 34) |
                         (sub1=='C' & str_sub(sub2, 1, 2) == 34), 1, 0)) %>% 
  
  # 유방암
  mutate(bc_p = ifelse((main1=='C' & str_sub(main2, 1, 2) == 50) |
                         (sub1=='C' & str_sub(sub2, 1, 2) == 50), 1, 0)) %>% 
  
  # T2DM
  mutate(t2dm_p = ifelse((main1=='E' & str_sub(main2, 1, 2) %in% 11:14) |
                           (sub1=='E' & str_sub(sub2, 1, 2) %in% 11:14), 1, 0)) %>% 
  
  # 고혈압
  mutate(hypertension_p = ifelse((main1=='I' & str_sub(main2, 1, 2) %in% 10:15) |
                                   (sub1=='I' & str_sub(sub2, 1, 2) %in% 10:15), 1, 0)) %>% 
  
  # 고지혈증
  mutate(dyslipidemia_p = ifelse((main1=='E' & str_sub(main2, 1, 2) == 78) |
                                   (sub1=='E' & str_sub(sub2, 1, 2) == 78), 1, 0)) %>% 
  
  # 천식
  mutate(asthma_p = ifelse((main1=='J' & str_sub(main2, 1, 2) %in% 45:46) |
                             (sub1=='J' & str_sub(sub2, 1, 2) %in% 45:46), 1, 0)) %>% 
  
  # 류머티스 관절염
  mutate(ra_p = ifelse((main1=='M' & str_sub(main2, 1, 2) == '05') |
                         (sub1=='M' & str_sub(sub2, 1, 2) == '05'), 1, 0)) %>% 
  
  # 만성폐쇄성폐질환
  mutate(copd_p = ifelse((main1=='J' & str_sub(main2, 1, 2) %in% 42:44) |
                           (sub1=='J' & str_sub(sub2, 1, 2) %in% 42:44), 1, 0)) %>% 
  
  # 결측 대치
  mutate_at(vars(stroke_p:copd_p), ~replace_na(., 0))



## 과거 기질환 추적기간 질환 처방약 투약여부
# 관찰년도 60T 로드
m160_df <- tibble()
for (y in base_year:base_year){
  path <- paste0(file_path, 'm160_', y, '.dta')
  tmp <- read_dta(path) %>% 
    rename_all(tolower) %>% 
    select(key_seq, gnl_nm_cd)
  
  m160_df <- bind_rows(m160_df, tmp)
  
  rm(tmp)
}






# m160에 약제를 처방 받았는 지 여부를 나타내는 더미변수 생성
m160_df2 <- m160_df %>% 
  mutate(p_dm_past = ifelse(gnl_nm_cd %in% dm_drug, 1, 0),
         p_hprts_past = ifelse(gnl_nm_cd %in% hprts_drug, 1, 0),
         p_asthma_past = ifelse(gnl_nm_cd %in% asthma_drug, 1, 0),
         p_copd_past = ifelse(gnl_nm_cd %in% copd_drug, 1, 0),
         p_ra_past = ifelse(gnl_nm_cd %in% ra_drug, 1, 0),
         p_statin_past = ifelse(gnl_nm_cd %in% statin_drug, 1, 0)) %>% 
  group_by(key_seq) %>% 
  summarise_at(vars(p_dm_past:p_statin_past), max)







# m120과 m160을 조인하여 관찰년도 중 개인이 한번이라도 약제를 처방받았는 지 여부를 나타내도록 groupby->max
m120_df3 <- m120_df2 %>% 
  left_join(m160_df2, by='key_seq') %>% 
  mutate_at(vars(contains('_past')), ~replace_na(., 0)) %>% 
  group_by(person_id) %>% 
  summarise_at(vars(stroke_p:p_statin_past), max)



write_csv(m120_df3, paste0(save_path, '3.past_disease_medicine_', base_year, '.csv'))


rm(m120)
rm(m120_df)
rm(m120_df2)
rm(m160_df)
rm(m160_df2)







##### 2. 타겟변수 전처리


# 약제코드 로드
dm_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\anti_dm_all.csv') %>% 
  rename(gnl_nm_cd = `Anti-diabeticmedicaiton(all)`) %>% 
  .$gnl_nm_cd
hprts_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\anti_hprts.csv') %>% 
  rename(gnl_nm_cd = `Anti-hypertensiveDrugs`) %>% 
  .$gnl_nm_cd
asthma_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\asthma.csv') %>% 
  .$gnl_nm_cd
copd_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\copd.csv') %>% 
  .$gnl_nm_cd
ra_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\ra.csv') %>% 
  .$gnl_nm_cd
statin_drug <- read_csv('C:\\Users\\wellxecon-r5\\Dropbox\\서상우\\drug_code\\statin.csv') %>% 
  rename(gnl_nm_cd = `Statin Drugs`) %>% 
  .$gnl_nm_cd




# 
# #### 2.1. VERSION1: 청구번호별 상병 판단
# target_ver1 <- tibble()
# for (y in (base_year+1):(base_year+track_windows)){
#   t0 <- Sys.time()
#   
#   # 연도별 m120 로드 후 전처리에 필요한 변수 추출
#   path_120 <- paste0(file_path, 'm120_', y, '.dta')
#   
#   m120 <- read_dta(path_120) %>% 
#     rename_all(tolower) %>% 
#     filter(person_id %in% eid) %>% 
#     
#     mutate(year = as.double(str_sub(recu_fr_dt, 1, 4))) %>% 
#     mutate(sub_sick = na_if(sub_sick, '')) %>% 
#     mutate(main1 = str_sub(main_sick, 1, 1),
#            main2 = str_sub(main_sick, 2, -1),
#            sub1 = str_sub(sub_sick, 1, 1),
#            sub2 = str_sub(sub_sick, 2, -1)) %>% 
#     select(person_id, key_seq, year, main1, main2, sub1, sub2, form_cd, vscn)
#   
#   
#   # 연도별 m160 로드 후 약제 처방 여부 변수 생성
#   path_160 <- paste0(file_path, 'm160_', y, '.dta')
#   m160 <- read_dta(path_160) %>% 
#     rename_all(tolower) %>% 
#     select(key_seq, gnl_nm_cd) %>% 
#     mutate(dm_prescrpt = ifelse(gnl_nm_cd %in% dm_drug, 1, 0),
#            hprts_prescrpt = ifelse(gnl_nm_cd %in% hprts_drug, 1, 0),
#            asthma_prescrpt = ifelse(gnl_nm_cd %in% asthma_drug, 1, 0),
#            copd_prescrpt = ifelse(gnl_nm_cd %in% copd_drug, 1, 0),
#            ra_prescrpt = ifelse(gnl_nm_cd %in% ra_drug, 1, 0),
#            statin_prescrpt = ifelse(gnl_nm_cd %in% statin_drug, 1, 0))
#   
#   
#   
#   # 청구번호별 상병, 입내원, 처방약 등을 활용하여 질병 발병 여부 타겟변수 생성
#   tmp <- m120 %>% 
#     left_join(m160, by='key_seq') %>% 
#     select(-gnl_nm_cd) %>% 
#     mutate_at(vars(dm_prescrpt:statin_prescrpt), ~replace_na(., 0)) %>% 
#     
#     
#     
#     # 뇌졸중
#     mutate(stroke_f = ifelse(((main1=='I' & str_sub(main2, 1, 2) %in% 60:69) |
#                                 (sub1=='I' & str_sub(sub2, 1, 2) %in% 60:69)) &
#                                (form_cd %in% c(2, 7, 10) & vscn >= 2), 1, 0)) %>% 
#     
#     # 심근경색
#     mutate(mi_f = ifelse(((main1=='I' & str_sub(main2, 1, 2) %in% 21:24) |
#                             (sub1=='I' & str_sub(sub2, 1, 2) %in% 21:24)) &
#                            (form_cd %in% c(2, 7, 10) & vscn >= 2), 1, 0)) %>% 
#     
#     # 골절
#     mutate(fracture_f = ifelse(((main1=='S' & str_sub(main2, 1, 2) %in% c(22, 32, 72)) | (main1=='M' & str_sub(main2, 1, 2) == 84) |
#                                   (sub1=='S' & str_sub(sub2, 1, 2) %in% c(22, 32, 72)) | (sub1=='M' & str_sub(sub2, 1, 2) == 84)) &
#                                  form_cd %in% c(2, 3, 7, 8, 10, 11), 1, 0)) %>% 
#     
#     # 치매
#     mutate(dementia_f = ifelse(((main1=='F' & str_sub(main2, 1, 2) %in% c('00', '01', '02', '03')) | (main1=='G' & str_sub(main2, 1, 2) %in% 30:32) |
#                                   (sub1=='F' & str_sub(sub2, 1, 2) %in% c('00', '01', '02', '03')) | (sub1=='G' & str_sub(sub2, 1, 2) %in% 30:32)) &
#                                  ((form_cd %in% c(2, 7, 10) & vscn >= 1) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 혈관성 치매
#     mutate(v_dementia_f = ifelse(((main1=='F' & str_sub(main2, 1, 2) %in% '01') | (sub1=='F' & str_sub(sub2, 1, 2) %in% '01')) &
#                                    ((form_cd %in% c(2, 7, 10) & vscn >= 1) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 우울증
#     mutate(depression_f = ifelse(((main1=='F' & str_sub(main2, 1, 2) %in% 32:33) |
#                                     (sub1=='F' & str_sub(sub2, 1, 2) %in% 32:33)) &
#                                    ((form_cd %in% c(2, 7, 10) & vscn >= 1) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 심부전
#     mutate(chf_f = ifelse(((main1=='I' & str_sub(main2, 1, 2) == 50) |
#                              (sub1=='I' & str_sub(sub2, 1, 2) == 50)) &
#                             form_cd %in% c(2, 7, 10), 1, 0)) %>% 
#     
#     # 폐렴
#     mutate(pneumonia_f = ifelse(((main1=='J' & str_sub(main2, 1, 2) %in% 12:17) |
#                                    (sub1=='J' & str_sub(sub2, 1, 2) %in% 12:17)) &
#                                   form_cd %in% c(2, 7, 10), 1, 0)) %>% 
#     
#     # 만성신장병
#     mutate(ckd_f = ifelse(((main1=='N' & str_sub(main2, 1, 3) %in% 183:185) |
#                              (sub1=='N' & str_sub(sub2, 1, 3) %in% 183:185)) &
#                             ((form_cd %in% c(2, 7, 10) & vscn >= 1) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 위암
#     mutate(gc_f = ifelse(((main1=='C' & str_sub(main2, 1, 2) == 16) |
#                             (sub1=='C' & str_sub(sub2, 1, 2) == 16)) &
#                            ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>%   
#     
#     # 대장암
#     mutate(cc_f = ifelse(((main1=='C' & str_sub(main2, 1, 2) %in% 18:21) |
#                             (sub1=='C' & str_sub(sub2, 1, 2) %in% 18:21)) &
#                            ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 간암
#     mutate(hc_f = ifelse(((main1=='C' & str_sub(main2, 1, 2) == 22) |
#                             (sub1=='C' & str_sub(sub2, 1, 2) == 22)) &
#                            ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 폐암
#     mutate(lc_f = ifelse(((main1=='C' & str_sub(main2, 1, 2) == 34) |
#                             (sub1=='C' & str_sub(sub2, 1, 2) == 34)) &
#                            ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 유방암
#     mutate(bc_f = ifelse(((main1=='C' & str_sub(main2, 1, 2) == 50) |
#                             (sub1=='C' & str_sub(sub2, 1, 2) == 50)) &
#                            ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # T2DM
#     mutate(t2dm_f = ifelse(((main1=='E' & str_sub(main2, 1, 2) %in% 11:14) |
#                               (sub1=='E' & str_sub(sub2, 1, 2) %in% 11:14)) &
#                              dm_prescrpt == 1, 1, 0)) %>% 
#     
#     # 고혈압
#     mutate(hypertension_f = ifelse(((main1=='I' & str_sub(main2, 1, 2) %in% 10:15) |
#                                       (sub1=='I' & str_sub(sub2, 1, 2) %in% 10:15)) &
#                                      hprts_prescrpt == 1, 1, 0)) %>% 
#     
#     # 고지혈증
#     mutate(dyslipidemia_f = ifelse(((main1=='E' & str_sub(main2, 1, 2) == 78) |
#                                       (sub1=='E' & str_sub(sub2, 1, 2) == 78)) &
#                                      statin_prescrpt == 1, 1, 0)) %>% 
#     
#     # 천식
#     mutate(asthma_f = ifelse(((main1=='J' & str_sub(main2, 1, 2) %in% 45:46) |
#                                 (sub1=='J' & str_sub(sub2, 1, 2) %in% 45:46)) &
#                                asthma_prescrpt == 1, 1, 0)) %>% 
#     
#     # 류머티스 관절염
#     mutate(ra_f = ifelse(((main1=='M' & str_sub(main2, 1, 2) == '05') |
#                             (sub1=='M' & str_sub(sub2, 1, 2) == '05')) &
#                            ra_prescrpt == 1, 1, 0)) %>% 
#     
#     # 만성폐쇄성폐질환
#     mutate(copd_f = ifelse(((main1=='J' & str_sub(main2, 1, 2) %in% 42:44) |
#                               (sub1=='J' & str_sub(sub2, 1, 2) %in% 42:44)) &
#                              copd_prescrpt == 1, 1, 0)) %>% 
#     
#     # 전체 암
#     mutate(cancer_f = ifelse(((main1=='C' & as.double(str_sub(main2, 1, 2)) %in% 0:99) |
#                                 (sub1=='C' & as.double(str_sub(sub2, 1, 2)) %in% 0:99)) &
#                                ((form_cd %in% c(2, 7, 10) & vscn >= 2) | (form_cd %in% c(3, 8, 11) & vscn >= 3)), 1, 0)) %>% 
#     
#     # 전체 뇌심혈관
#     mutate(vessel_f = ifelse(((main1=='I' & as.double(str_sub(main2, 1, 2)) %in% c(20:25, 60:69)) |
#                                 (sub1=='I' & as.double(str_sub(sub2, 1, 2)) %in% c(20:25, 60:69))) &
#                                (form_cd %in% c(2, 7, 10) & vscn >= 2), 1, 0)) %>% 
#     
#     
#     # 결측 대치
#     mutate_at(vars(stroke_f:vessel_f), ~replace_na(., 0))
#   
#   
#   # ID별 groupby -> max
#   tmp <- tmp %>% 
#     group_by(person_id, year) %>% 
#     summarise_at(vars(stroke_f:vessel_f), ~max(.)) %>% 
#     ungroup()
#   
#   target_ver1 <- bind_rows(target_ver1, tmp)
#   
#   rm(m120)
#   rm(m160)
#   rm(tmp)
#   
#   t1 <- Sys.time()
#   cat('Running time:', t1-t0, 'm\n', sep='')
# }
# 
# 
# 
# 
# 
# write_csv(target_ver1, paste0(save_path, '4_1.target_tb_v1_', base_year, '.csv'), na="")
# 
# 




#### 2.2. VERSION2: 연단위 상병 판단
# 연도별 m120 로드 후 전처리에 필요한 변수 추출
target_ver2 <- tibble()
for(y in (base_year+1):(base_year+track_windows)){
  t0 <- Sys.time()
  
  # m120 
  path_120 <- paste0(file_path, 'm120_', y, '.dta')
  m120 <- read_dta(path_120) %>% 
    rename_all(tolower)
  
  m120_2 <- p_final %>% 
    filter(person_id %in% eid) %>% 
    select(person_id) %>% 
    left_join(m120, by='person_id') %>% 
    
    mutate(year = as.double(str_sub(recu_fr_dt, 1, 4))) %>% 
    mutate(sub_sick = na_if(sub_sick, '')) %>% 
    select(person_id, key_seq, year, main_sick, sub_sick, ykiho_id, recu_fr_dt, form_cd, vscn)
  
  
  
  # 연도별 m160 로드 후 약제 처방 여부 변수 생성
  path_160 <- paste0(file_path, 'm160_', y, '.dta')
  m160 <- read_dta(path_160) %>% 
    rename_all(tolower) %>% 
    select(key_seq, gnl_nm_cd) %>% 
    mutate(dm_prescrpt = ifelse(gnl_nm_cd %in% dm_drug, 1, 0),
           hprts_prescrpt = ifelse(gnl_nm_cd %in% hprts_drug, 1, 0),
           asthma_prescrpt = ifelse(gnl_nm_cd %in% asthma_drug, 1, 0),
           copd_prescrpt = ifelse(gnl_nm_cd %in% copd_drug, 1, 0),
           ra_prescrpt = ifelse(gnl_nm_cd %in% ra_drug, 1, 0),
           statin_prescrpt = ifelse(gnl_nm_cd %in% statin_drug, 1, 0))
  
  
  
  tmp <- m120_2 %>% 
    left_join(m160, by='key_seq') %>% 
    
    
    mutate(year = replace_na(year, y)) %>%
    mutate_at(vars(dm_prescrpt:statin_prescrpt), ~replace_na(., 0)) %>% 
    
    group_by(person_id, year, main_sick, sub_sick, ykiho_id, recu_fr_dt, form_cd) %>% 
    summarise(vscn = max(vscn),
              dm_prescrpt = max(dm_prescrpt),
              hprts_prescrpt = max(hprts_prescrpt),
              asthma_prescrpt = max(asthma_prescrpt),
              copd_prescrpt = max(copd_prescrpt),
              ra_prescrpt = max(ra_prescrpt),
              statin_prescrpt = max(statin_prescrpt)) %>% 
    ungroup() %>% 
    
    mutate(main_sick = replace_na(main_sick, 'Na')) %>% 
    mutate(sub_sick = ifelse(str_sub(main_sick, 1, 3) == str_sub(sub_sick, 1, 3), NA_character_, sub_sick))
  
  
  
  
  
  tmp2 <- tmp %>% 
    gather('tmp', 'sick_code', main_sick, sub_sick) %>%
    select(-tmp) %>% 
    filter(!is.na(sick_code)) %>% 
    mutate(sick_code2 = str_sub(sick_code, 1, 3),
           in_day = ifelse(form_cd %in% c(2, 7, 10), vscn, 0),
           out_day = ifelse(form_cd %in% c(3, 8, 11), vscn, 0),
           form_cd = case_when(form_cd %in% c(2, 7, 10) ~ 'in',
                               form_cd %in% c(3, 8, 11) ~ 'out',
                               TRUE ~ NA_character_))
  
  
  
  
  prescrpt <- tmp2 %>% 
    group_by(person_id, year) %>% 
    summarise_at(vars(ends_with('_prescrpt')), max) %>% 
    ungroup()
  
  
  
  tmp3 <- tmp2 %>% 
    group_by(person_id, year, form_cd, sick_code2) %>% 
    summarise(in_day = sum(in_day),
              out_day = sum(out_day)) %>% 
    ungroup() %>% 
    left_join(prescrpt, by=c('person_id', 'year')) %>% 
    
    mutate(sick1 = str_sub(sick_code2, 1, 1),
           sick2 = str_sub(sick_code2, 2, -1))
  
  
  
  
  
  
  tmp4 <- tmp3 %>% 
    
    # 뇌졸중 1
    mutate(stroke_f = ifelse((sick1=='I' & as.double(str_sub(sick2, 1, 2)) %in% 60:69) &
                               (in_day >= 2), 1, 0)) %>% 
    
    # 심근경색 2
    mutate(mi_f = ifelse((sick1=='I' & as.double(str_sub(sick2, 1, 2)) %in% 21:24) &
                           (in_day >= 2), 1, 0)) %>% 
    
    # 골절 3 
    mutate(fracture_f = ifelse(((sick1=='S' & as.double(str_sub(sick2, 1, 2)) %in% c(22, 32, 72)) | (sick1=='M' & as.double(str_sub(sick2, 1, 2)) == 84)) &
                                 ((form_cd=='in' & in_day>=1) | (form_cd=='out' & out_day>=1)), 1, 0)) %>% 
    
    # 치매 4
    mutate(dementia_f = ifelse(((sick1=='F' & as.double(str_sub(sick2, 1, 2)) %in% 0:3) | (sick1=='G' & as.double(str_sub(sick2, 1, 2)) %in% 30:32)) & 
                                 ((in_day >= 1) | (out_day >=3)), 1, 0)) %>% 
    
    # 혈관성 치매 5
    mutate(v_dementia_f = ifelse((sick1=='F' & str_sub(sick2, 1, 2) %in% '01') &
                                   ((in_day >= 1) | (out_day >= 3)), 1, 0)) %>% 
    
    # 우울증 6
    mutate(depression_f = ifelse((sick1=='F' & as.double(str_sub(sick2, 1, 2)) %in% 32:33) &
                                   ((in_day >= 1) | (out_day >= 3)), 1, 0)) %>% 
    
    # 심부전 7
    mutate(chf_f = ifelse((sick1=='I' & as.double(str_sub(sick2, 1, 2)) == 50) &
                            (form_cd == 'in' & in_day >= 1), 1, 0)) %>% 
    
    # 폐렴 8
    mutate(pneumonia_f = ifelse((sick1=='J' & as.double(str_sub(sick2, 1, 2)) %in% 12:17) &
                                  (form_cd == 'in' & in_day >= 1), 1, 0)) %>% 
    
    # 만성신장병 9
    mutate(ckd_f = ifelse((sick1=='N' & as.double(str_sub(sick2, 1, 2)) == 18) &
                            ((in_day >= 1) | (out_day >= 3)), 1, 0)) %>% 
    
    # 위암 10
    mutate(gc_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) == 16) &
                           ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>%   
    
    # 대장암 11
    mutate(cc_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) %in% 18:21) &
                           ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>% 
    
    # 간암 12
    mutate(hc_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) == 22) &
                           ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>% 
    
    # 폐암 13
    mutate(lc_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) == 34) &
                           ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>% 
    
    # 유방암 14
    mutate(bc_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) == 50) &
                           ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>% 
    
    # T2DM 15
    mutate(t2dm_f = ifelse((sick1=='E' & as.double(str_sub(sick2, 1, 2)) %in% 11:14) &
                             dm_prescrpt == 1, 1, 0)) %>% 
    
    # 고혈압 16
    mutate(hypertension_f = ifelse((sick1=='I' & as.double(str_sub(sick2, 1, 2)) %in% 10:15) &
                                     hprts_prescrpt == 1, 1, 0)) %>% 
    
    # 고지혈증 17
    mutate(dyslipidemia_f = ifelse((sick1=='E' & as.double(str_sub(sick2, 1, 2)) == 78) &
                                     statin_prescrpt == 1, 1, 0)) %>% 
    
    # 천식 18
    mutate(asthma_f = ifelse((sick1=='J' & as.double(str_sub(sick2, 1, 2)) %in% 45:46) &
                               asthma_prescrpt == 1, 1, 0)) %>% 
    
    # 류머티스 관절염 19
    mutate(ra_f = ifelse((sick1=='M' & as.double(str_sub(sick2, 1, 2)) == 5) &
                           ra_prescrpt == 1, 1, 0)) %>% 
    
    # 만성폐쇄성폐질환 20
    mutate(copd_f = ifelse((sick1=='J' & as.double(str_sub(sick2, 1, 2)) %in% 42:44) &
                             copd_prescrpt == 1, 1, 0)) %>% 
    
    
    # 전체 암
    mutate(cancer_f = ifelse((sick1=='C' & as.double(str_sub(sick2, 1, 2)) %in% 0:99) &
                               ((in_day >= 2) | (out_day >= 3)), 1, 0)) %>% 
    
    # 전체 뇌심혈관
    mutate(vessel_f = ifelse((sick1=='I' & as.double(str_sub(sick2, 1, 2)) %in% c(20:25, 60:69)) &
                               (in_day >= 2), 1, 0)) %>% 
    
    
    # 결측 대치
    mutate_at(vars(stroke_f:vessel_f), ~replace_na(., 0))
  
  
  
  
  death <- p_final %>% 
    filter(dth_year == y) %>% 
    select(person_id, dth_year) %>% 
    mutate(dth_year = 1)
  
  
  tmp5 <- tmp4 %>% 
    group_by(person_id, year) %>% 
    summarise_at(vars(stroke_f:vessel_f), max) %>% 
    ungroup() %>% 
    
    left_join(death, by='person_id') %>% 
    mutate(dth_year = replace_na(dth_year, 0)) %>% 
    
    mutate(mvp1 = apply(select(., cancer_f:dth_year, t2dm_f), 1, function(x) max(x)))
  
  
  
  target_ver2 <- bind_rows(target_ver2, tmp5)
  
  t1 <- Sys.time()
  cat('Running time: ', t1-t0, 'm\n', sep='')
  
  
  rm(tmp)
  rm(tmp2)
  rm(tmp3)
  rm(tmp4)
  rm(tmp5)
  rm(m120)
  rm(m120_2)
  rm(m160)
  rm(prescrpt)
  
}

target_ver2_final <- target_ver2 %>% 
  group_by(person_id) %>% 
  summarise_at(vars(stroke_f:mvp1), max) %>% 
  ungroup()


write_csv(target_ver2_final, paste0(save_path, '4_2.target_tb_v2_', base_year, '.csv'), na="")



