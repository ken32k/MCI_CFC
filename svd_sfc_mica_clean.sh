#!/bin/bash
###
 # @Author: Ken32g ken32k@163.com
 # @Date: 2024-07-01 14:01:07
 # @LastEditors: Ken32g ken32k@163.com
 # @LastEditTime: 2024-07-02 14:19:38
 # @FilePath: /csvd-sfc/svd_sfc_mica_clean.sh
 # @Description: Clean MICApipe output files
### 

sub=$1
bids=${CSVD_DIR}/csvd_mica_bids
out=${CSVD_DIR}/csvd_mica_out
micapipe_cleanup -bids $bids -out $out -sub $sub -proc_dwi -SC