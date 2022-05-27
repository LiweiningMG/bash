#!/usr/bin/bash
## 根据提供的标记物理位置进行LD连锁分析
## 需要提供的信息：
##   基因组区域文件，前四列为Trait CHR S_POS E_POS，其他列随意，有标题行

while getopts "hI:B:P:O:D:F:" optname
do
    case "$optname" in
      "I")
        snpf=${OPTARG}
        ;;
      "B")
        bfile=${OPTARG}
        ;;
      "P")
        pre=${OPTARG}
        ;;
      "O")
        out=${OPTARG}
        ;;
      "D")
        dir=${OPTARG}
        ;;
      "h")
        echo "Usage: ./haploview_run.sh -I snp_pos.txt -B bfile_prefix -P pheno_label -D directory -O prefix"
        exit -1
        ;;
      ":")
        echo "No argument value for option $OPTARG"
        ;;
      "?")
        echo "Unknown option $OPTARG"
        ;;
      *)
        echo "Unknown error while processing options"
        ;;
    esac
    #echo "option index is $OPTIND"
done

## 软件加载
module load PLINK/1.90

## 文件/文件夹
pedi=/storage1/active/liwn/tempWork/FCR_GWAS/ped_BGYD.txt

## 软件/脚本
Haploview=/storage1/active/liwn/Mysoftware/Haploview.jar
haploview_format=/storage1/active/liwn/code/R/LD/haploview_file_format.R
snp_up_down=/storage1/active/liwn/code/R/LD/snp_up_downstream.R
block_check=/storage1/active/liwn/code/R/LD/block_check.R

## 标记位点文件
[[ ! -f ${snpf} ]] && echo "snp file ${snpf} not found" && exit -1

rows=`wc -l < ${snpf}`
for r in `seq 2 ${rows}`; do
  trait=`sed -n ${r}p ${snpf} | awk '{print $1}'`
  chr=`sed -n ${r}p ${snpf} | awk '{print $2}'`
  start_pos=`sed -n ${r}p ${snpf} | awk '{print $3}'`
  end_pos=`sed -n ${r}p ${snpf} | awk '{print $4}'`
  snp_id=`sed -n ${r}p ${snpf} | awk '{print $5}'`
  nSNP=`sed -n ${r}p ${snpf} | awk '{print $6}'`
  echo ${trait} ${chr} ${start_pos} ${end_pos} ${nSNP} ${snp_id}
  
  # ## 单个snp的QTL跳过
  # [[ ${nSNP} -le 1 ]] && echo "$((r -1))th QTL contains only 1 SNP, skip..." && continue
  
  # ## 按qtl区间提取位点，同时将碱基转为1=A; 2=C; 3=G; 4=T; 缺失碱基用0表示
  # plink --bfile ${bfile} \
  #   --chr ${chr} \
  #   --from-bp ${start_pos} --to-bp ${end_pos} \
  #   --allele1234 \
  #   --recode --out ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp
  
  # ## 只提取显著标记
  sed -n ${r}p ${snpf} | awk '{print $5}' | tr ',' '\n' > SNPs_extract.txt
  
  ## 提取每个qtl中的显著标记及上下游10个snp（在每个qtl上下游加了10个SNP为了防止tagSNPs边上的连锁SNP丢)
  $snp_up_down \
    --file ${bfile} \
    --sig ${snp_id} \
    --near 10 \
    --sigRank ${dir}/sig_index.txt \
    --out ${dir}/SNPs_extract.txt
  
  plink --bfile ${bfile} \
      --extract ${dir}/SNPs_extract.txt \
      --allele1234 \
      --recode --out ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp && \
      rm ${dir}/SNPs_extract.txt
  
  ## 转换为haploview所需格式
  $haploview_format \
    --file ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp \
    --pedi ${pedi} \
    --out ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp
  
  ## haploview作图
  java -jar $Haploview \
    -skipcheck -n \
    -pedfile ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp.ped \
    -info ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp.info \
    -dprime \
    -svg \
    -nogui \
    -quiet \
    -blockoutput GAB \
    -out ${dir}/${pre}
  
  ## 检查block中是否有显著标记
  hit=F
  for i in `cat ${dir}/sig_index.txt`; do
    block=`grep "MARKERS" ${dir}/${pre}.GABRIELblocks | grep ${i} | wc -l`
    if [[ ${block} -gt 0 ]]; then
      hit=T
    fi
  done
  
  if [[ ${hit} == 'T' ]]; then
    ## 更换图片背景颜色
    sed -i 's/rgb(212,208,200)/none/g' ${dir}/${pre}.LD.SVG
    
    ## 改名(haploview有bug)
    mv ${dir}/${pre}.LD.SVG ${dir}/${pre}_${trait}_${chr}_${start_pos}_${nSNP}snp.SVG
    mv ${dir}/${pre}.LD ${dir}/${pre}_${trait}_${chr}_${start_pos}_${nSNP}snp.LD
    
    echo LD block present in qtl $((r - 1))!    
  else
    rm ${pre}.LD.SVG ${pre}.LD
    echo No LD block in qtl $((r - 1))!
  fi
  
  ## 删除中间文件
  
  rm  ${dir}/${pre}_${trait}_${chr}_${start_pos}_${end_pos}_${nSNP}snp.*
done

[[ -f ${dir}/${pre}.GABRIELblocks ]] && rm  ${dir}/${pre}.GABRIELblocks
[[ -f ${dir}/sig_index.txt ]] && rm  ${dir}/sig_index.txt

## debug
snpf=${root}/GEMMA/${b}/${p}/${b}_${p}_region.txt
bfile=${root}/chip/${b}_beagle_11.1_pca_auto_qc
pre=${b}_${p}
dir=${root}/haploview/${b}/${p}
