import pandas as pd
import numpy as np
import argparse

def calculate_slope(y):
    """计算给定窗口的斜率"""
    x = np.array(range(0, y.shape[0]))  # x坐标从0到窗口大小-1
    slope, intercept = np.polyfit(x, y, 1)
    return slope

def load_data(gene_bed_path, tu_file_path, coverage_file_path):
    """加载输入文件"""
    gene_bed = pd.read_csv(gene_bed_path, sep='\t', header=None)
    tu_df = pd.read_csv(tu_file_path)
    coverage_df = pd.read_csv(coverage_file_path, sep='\t', header=None)
    return gene_bed, tu_df, coverage_df

def extract_gene_info(tu_df, gene_bed):
    """从TU数据和BED文件中提取基因信息"""
    result_list = []
    for idx, row in tu_df.iterrows():
        gene_list = row.iloc[0].split('|')
        strand = None
        for gene in gene_list:
            gene_info = gene_bed[gene_bed[3] == gene]
            if not gene_info.empty:  # 确保gene_info不为空
                if strand is None:
                    strand = gene_info.iloc[0][5]
                elif strand != gene_info.iloc[0][5]:
                    continue
        gene_df = gene_bed[gene_bed[3].isin(gene_list)]
        if not gene_df.empty:  # 确保gene_df不为空
            CDS_start = gene_df[1].min()
            CDS_end = gene_df[2].max()
            result_list.append([gene_df.iloc[0][0], CDS_start, CDS_end,row.iloc[0], strand])
    return result_list

def calculate_utr_regions(result_list, coverage_df, tss_range=150, window_size=6):
    """计算UTR起点和终点"""
    result_df_list = []
    for item in result_list:
        CDS_start = [item[1] - tss_range, item[1]]
        CDS_end = [item[2], item[2] + tss_range]

        # UTR start
        start_df = coverage_df[(coverage_df[1] >= CDS_start[0]) & (coverage_df[1] <= CDS_start[1])&(coverage_df[0]==item[0])]
        if not start_df.empty:  # 确保窗口数据足够
            slope = start_df[2].rolling(window=window_size).apply(calculate_slope)
            slope = slope[slope>0]
            if slope.empty:
                UTR_start = None
            else:
                UTR_start = slope.idxmax()
        else:
            UTR_start = None  # 如果数据不足，返回None

        # UTR end
        end_df = coverage_df[(coverage_df[1] >= CDS_end[0]) & (coverage_df[1] <= CDS_end[1])&(coverage_df[0]==item[0])]
        if not end_df.empty:   # 确保窗口数据足够
            slope = end_df[2].rolling(window=window_size).apply(calculate_slope)
            slope = slope[slope < 0]
            if slope.empty:
                UTR_end = None
            else:
                UTR_end = slope.idxmin()
        else:
            UTR_end = None  # 如果数据不足，返回None

        item.extend([UTR_start, UTR_end])
        result_df_list.append(item)
    return result_df_list

def save_results(result_df_list, output_path):
    """保存结果到CSV文件"""
    result_df = pd.DataFrame(result_df_list)
    result_df.columns = ['Chromosome', 'CDS start', 'CDS end','TU', 'strand', 'UTR start', 'UTR end']
    result_df.to_csv(output_path, index=False)

def main():
    """主函数：处理命令行参数并协调各个模块"""
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Calculate TSS and UTR regions from gene and coverage data.")
    parser.add_argument('--gene_bed', type=str, default="/t1/zhguo/Data/dx_RNA/ref/Ecoli_gene.bed", help="Path to the gene BED file")
    parser.add_argument('--tu_file', type=str, default='/t1/zhguo/Data/dx_RNA/k12_od02//longest_tu.csv', help="Path to the TU CSV file")
    parser.add_argument('--coverage_file', type=str, default='/t1/zhguo/Data/dx_RNA/k12_od02//coverage.txt', help="Path to the coverage file")
    parser.add_argument('--tss_range', type=int, default=150,help="Region of distance from CDS where to find TSS regions")
    parser.add_argument('--slope_windows', type=int, default=6,help="Windows size of slope calculation")
    parser.add_argument('--output', type=str, default='/t1/zhguo/Data/dx_RNA/k12_od02/TSS.csv', help="Path to the output CSV file")
    args = parser.parse_args()

    # 加载数据
    gene_bed, tu_df, coverage_df = load_data(args.gene_bed, args.tu_file, args.coverage_file)

    # 提取基因信息
    result_list = extract_gene_info(tu_df, gene_bed)

    # 计算UTR区域
    result_df_list = calculate_utr_regions(result_list, coverage_df, args.tss_range, args.slope_windows)

    # 保存结果
    save_results(result_df_list, args.output)

    print("Processing completed. Output saved to:", args.output)

if __name__ == "__main__":
    main()
