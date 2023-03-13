#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: drawStressed.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
'''

import csv
import os
import re
import time

os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2, 40).__str__()
import pickle
import random
import warnings

import anndata as ad
import cv2
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import TESLA as tesla
from IPython.display import Image
from scipy import stats
from scipy.sparse import issparse


def draw(matrixH5Path, spCSVPath, imgPath, tag, res, outDir):

    # matrixH5Path = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/TESLA/data/for_Kevin/GSE203612_RAW_2/GSM6177599/GSM6177599_NYU_BRCA0_Vis_processed_filtered_feature_bc_matrix.h5"
    # spCSVPath = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/TESLA/data/for_Kevin/GSE203612_RAW_2/GSM6177599/GSM6177599_NYU_BRCA0_Vis_processed_spatial_tissue_positions_list.csv"
    # imgPath = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/TESLA/data/for_Kevin/GSE203612_RAW_2/GSM6177599/GSM6177599_NYU_BRCA0_Vis_processed_spatial_tissue_hires_image.png"
    # tag = "GSE203612"
    # outDir = "/rsrch3/scratch/genomic_med/ychu2/projects/p1review/TESLA/result/1_draw_stressed_T_find_res_black/GSE203612_RAW_2/GSM6177599/50"
    # res = 50

    adata = sc.read_10x_h5(matrixH5Path)
    spatial = pd.read_csv(spCSVPath,
                          sep=",",
                          header=None,
                          na_filter=False,
                          index_col=0)
    img = cv2.imread(imgPath)

    adata.obs["x1"] = spatial[1]
    adata.obs["x2"] = spatial[2]
    adata.obs["x3"] = spatial[3]
    adata.obs["x4"] = spatial[4]
    adata.obs["x5"] = spatial[5]
    adata.obs["array_x"] = adata.obs["x2"]
    adata.obs["array_y"] = adata.obs["x3"]
    adata.obs["pixel_x"] = adata.obs["x4"]
    adata.obs["pixel_y"] = adata.obs["x5"]

    adata = adata[adata.obs["x1"] == 1]
    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")

    counts = adata
    resize_factor = 2000 / np.min(img.shape[0:2])
    resize_width = int(img.shape[1] * resize_factor)
    resize_height = int(img.shape[0] * resize_factor)
    counts.var.index = [i.upper() for i in counts.var.index]
    counts.var_names_make_unique()
    counts.raw = counts
    sc.pp.log1p(counts)  # impute on log scale
    if issparse(counts.X):
        counts.X = counts.X.A

    print("Step 1: Gene expression enhancement ")

    cnt = tesla.cv2_detect_contour(img, apertureSize=5, L2gradient=True)
    binary = np.zeros((img.shape[0:2]), dtype=np.uint8)
    cv2.drawContours(binary, [cnt], -1, (1), thickness=-1)

    enhanced_exp_adata = tesla.imputation(img=img,
                                          raw=counts,
                                          cnt=cnt,
                                          genes=counts.var.index.tolist(),
                                          shape="None",
                                          res=res,
                                          s=1,
                                          k=2,
                                          num_nbs=10)
    print("Imputation done")

    TGenes = ["CD3D", "CD3E", "CD3G"]

    genes = list(set([i for i in TGenes if i in enhanced_exp_adata.var.index]))
    pred_refined_T, target_clusters_T, c_m_T = tesla.annotation(
        img=img,
        binary=binary,
        sudo_adata=enhanced_exp_adata,
        genes=genes,
        resize_factor=resize_factor,
        num_required=2,
        target_size="small")
    topClassesT = [i for i, v in c_m_T if v == max([v2 for _, v2 in c_m_T])]

    ret_img = tesla.visualize_annotation(img=img,
                                         binary=binary,
                                         resize_factor=resize_factor,
                                         pred_refined=pred_refined_T,
                                         target_clusters=topClassesT,
                                         c_m=c_m_T)
    cv2.imwrite(outDir + '/CD3-' + tag + '.jpg', ret_img)
    Image(filename=outDir + '/CD3-' + tag + '.jpg')

    ret_img = tesla.visualize_annotation(img=img,
                                         binary=binary,
                                         resize_factor=resize_factor,
                                         pred_refined=pred_refined_T,
                                         target_clusters=target_clusters_T,
                                         c_m=c_m_T)
    cv2.imwrite(outDir + '/CD3-all-' + tag + '.jpg', ret_img)
    Image(filename=outDir + '/CD3-all-' + tag + '.jpg')

    genes = ["HSPA1A", "HSPA1B", "CD3D", "CD3E", "CD3G"]
    genes = list(set([i for i in genes if i in enhanced_exp_adata.var.index]))
    pred_refined_Ts, target_clusters_Ts, c_m_Ts = tesla.annotation(
        img=img,
        binary=binary,
        sudo_adata=enhanced_exp_adata,
        genes=genes,
        resize_factor=resize_factor,
        num_required=5,
        target_size="small")
    ret_img = tesla.visualize_annotation(img=img,
                                         binary=binary,
                                         resize_factor=resize_factor,
                                         pred_refined=pred_refined_Ts,
                                         target_clusters=target_clusters_Ts,
                                         c_m=c_m_Ts)
    cv2.imwrite(outDir + '/HSP-' + tag + '.jpg', ret_img)
    Image(filename=outDir + '/HSP-' + tag + '.jpg')

    CD3Position = np.isin(pred_refined_T, target_clusters_T)
    HSPPosition = np.isin(pred_refined_Ts, target_clusters_Ts)
    newCluster = max([i for i, j in c_m_Ts]) + 1
    c_m_Ts_filtered = c_m_Ts[0:len(c_m_Ts)] + [(newCluster, 0.0)]
    pred_refined_Ts_filtered = pred_refined_Ts
    pred_refined_Ts_filtered[np.logical_and(np.invert(CD3Position),
                                            HSPPosition)] = newCluster

    topClassesTs = [
        i for i, v in c_m_Ts_filtered
        if v == max([v2 for _, v2 in c_m_Ts_filtered])
    ]

    ret_img = tesla.visualize_annotation(img=img,
                                         binary=binary,
                                         resize_factor=resize_factor,
                                         pred_refined=pred_refined_Ts_filtered,
                                         target_clusters=topClassesTs,
                                         c_m=c_m_Ts_filtered)
    cv2.imwrite(outDir + '/CD3HSP-' + tag + '.jpg', ret_img)
    Image(filename=outDir + '/CD3HSP-' + tag + '.jpg')

    ret_img = tesla.visualize_annotation(img=img,
                                         binary=binary,
                                         resize_factor=resize_factor,
                                         pred_refined=pred_refined_Ts_filtered,
                                         target_clusters=target_clusters_Ts,
                                         c_m=c_m_Ts_filtered)
    cv2.imwrite(outDir + '/CD3HSP-all-' + tag + '.jpg', ret_img)
    Image(filename=outDir + '/CD3HSP-all-' + tag + '.jpg')

    hypoxia_genes = [
        "HIF1A", "EPAS1", "RELA", "RELB", "NFKB1", "NFKB2", "NFE2L2", "CREB1"
    ]
    hypoxia_genes = list(
        set([i for i in hypoxia_genes if i in enhanced_exp_adata.var.index]))

    for genei in hypoxia_genes:
        try:
            pred_refined_h, target_clusters_h, c_m_h = tesla.annotation(
                img=img,
                binary=binary,
                sudo_adata=enhanced_exp_adata,
                genes=[genei],
                resize_factor=resize_factor,
                num_required=1,
                target_size="small")
            topClasses_h = [
                i for i, v in c_m_h if v == max([v2 for _, v2 in c_m_h])
            ]

            ret_img, black_ret_img = visualize_annotation_double(
                pred_refined_1=pred_refined_Ts_filtered,
                target_clusters_1=target_clusters_Ts,
                c_m_1=c_m_Ts_filtered,
                pred_refined_2=pred_refined_h,
                target_clusters_2=topClasses_h,
                c_m_2=c_m_h,
                img=img,
                binary=binary,
                resize_factor=resize_factor)
            cv2.imwrite(outDir + '/hypoxia-' + genei + '.jpg', ret_img)
            Image(filename=outDir + '/hypoxia-' + genei + '.jpg')
            cv2.imwrite(outDir + '/hypoxia-black-' + genei + '.jpg',
                        black_ret_img)
            Image(filename=outDir + '/hypoxia-black-' + genei + '.jpg')

            ret_img, black_ret_img = visualize_annotation_double(
                pred_refined_1=pred_refined_Ts_filtered,
                target_clusters_1=target_clusters_Ts,
                c_m_1=c_m_Ts_filtered,
                pred_refined_2=pred_refined_h,
                target_clusters_2=target_clusters_h,
                c_m_2=c_m_h,
                img=img,
                binary=binary,
                resize_factor=resize_factor)
            cv2.imwrite(outDir + '/hypoxia-all-' + genei + '.jpg', ret_img)
            Image(filename=outDir + '/hypoxia-all-' + genei + '.jpg')
            cv2.imwrite(outDir + '/hypoxia-all-black-' + genei + '.jpg',
                        black_ret_img)
            Image(filename=outDir + '/hypoxia-all-black-' + genei + '.jpg')

        except:
            print(genei, " failed")


def visualize_annotation_double(
    pred_refined_1,
    target_clusters_1,
    c_m_1,
    pred_refined_2,
    target_clusters_2,
    c_m_2,
    img,
    binary,
    resize_factor,
    cnt_color_1=clr.LinearSegmentedColormap.from_list('red',
                                                      ["#EAE7CC", '#BA0000'],
                                                      N=256),
    cnt_color_2=clr.LinearSegmentedColormap.from_list('green',
                                                      ["#E7EACC", '#00BA00'],
                                                      N=256),
    cnt_color_3=clr.LinearSegmentedColormap.from_list('blue',
                                                      ["#CCEAE7", '#5555EE'],
                                                      N=256)):
    resize_width = int(img.shape[1] * resize_factor)
    resize_height = int(img.shape[0] * resize_factor)
    binary_resized = cv2.resize(binary, (resize_width, resize_height),
                                interpolation=cv2.INTER_AREA)
    background = cv2.resize(img, (resize_width, resize_height),
                            interpolation=cv2.INTER_AREA)
    white_background = np.ones(shape=np.shape(background)) * 255

    ret_img = (background.copy()).astype(np.uint8)
    white_ret_img = (white_background.copy()).astype(np.uint8)

    alpha = 0.8
    #Whiten
    white_ratio = 0.5
    ret_img = ret_img * (1 - white_ratio) + np.array([255, 255, 255
                                                      ]) * (white_ratio)

    target_img_1_total = (
        1 * (np.isin(pred_refined_1, target_clusters_1))).reshape(
            resize_height, resize_width)
    target_img_1_total[binary_resized == 0] = 0

    for i in range(len(target_clusters_1)):
        color_1 = ((np.array(cnt_color_1(int(
            c_m_1[i][1] / c_m_1[0][1] * 255)))[0:3]) * 255).astype(int)[::-1]

        target_img_1 = (1 * (pred_refined_1 == target_clusters_1[i])).reshape(
            resize_height, resize_width)

        target_img_1[binary_resized == 0] = 0

        ret_img[target_img_1 != 0] = ret_img[target_img_1 != 0] * (
            1 - alpha) + np.array(color_1) * (alpha)

        white_ret_img[target_img_1 != 0] = np.array(color_1)

    alpha2 = 0.5
    alpha3 = 0.2
    for i in range(len(target_clusters_2)):
        color_2 = ((np.array(cnt_color_2(int(
            c_m_2[i][1] / c_m_2[0][1] * 255)))[0:3]) * 255).astype(int)[::-1]
        color_3 = ((np.array(cnt_color_3(int(
            c_m_2[i][1] / c_m_2[0][1] * 255)))[0:3]) * 255).astype(int)[::-1]

        target_img_2 = (1 * (pred_refined_2 == target_clusters_2[i])).reshape(
            resize_height, resize_width)
        target_img_2[binary_resized == 0] = 0

        notCovered = np.logical_and(target_img_2 != 0, target_img_1_total == 0)

        ret_img[notCovered] = ret_img[notCovered] * (
            1 - alpha) + np.array(color_2) * (alpha)
        white_ret_img[notCovered] = np.array(color_3)

        covered = np.logical_and(target_img_2 != 0, target_img_1_total != 0)
        ret_img[covered] = ret_img[covered] * (
            1 - alpha2) + np.array(color_2) * (alpha2)
        white_ret_img[covered] = white_ret_img[covered] * (
            1 - alpha3) + np.array(color_3) * (alpha3)

    return ret_img, white_ret_img


def main():
    import argparse
    parser = argparse.ArgumentParser(description='draw stressed figure')
    parser.add_argument('--matrixH5Path',
                        dest='matrixH5Path',
                        help='matrixH5Path')
    parser.add_argument('--spCSVPath', dest='spCSVPath', help='spCSVPath')
    parser.add_argument('--imgPath', dest='imgPath', help='imgPath')
    parser.add_argument('--tag', dest='tag', help='tag')
    parser.add_argument('--res', dest='res', help='res')
    parser.add_argument('--outDir', dest='outDir', help='outDir')
    args = parser.parse_args()
    draw(args.matrixH5Path, args.spCSVPath, args.imgPath, args.tag,
         float(args.res), args.outDir)


if __name__ == '__main__':
    main()
