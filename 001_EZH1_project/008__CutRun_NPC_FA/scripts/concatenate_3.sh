#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


cat input_raw_Novogene/NPC_OE_KO_IGG_CKDL240001582-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_OE_KO_IGG_CKDL240001582-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_OE_KO_IGG_L14_2.fq.gz
cat input_raw_Novogene/NPC_OE_KO_SUZ12_CKDL240001588-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_OE_KO_SUZ12_CKDL240001588-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_OE_KO_SUZ12_L14_1.fq.gz
cat input_raw_Novogene/NPC_OE_KO_SUZ12_CKDL240001588-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_OE_KO_SUZ12_CKDL240001588-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_OE_KO_SUZ12_L14_2.fq.gz
cat input_raw_Novogene/NPC_WT_EZH1_CS_CKDL240001572-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_WT_EZH1_CS_CKDL240001572-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_WT_EZH1_CS_L14_1.fq.gz
cat input_raw_Novogene/NPC_WT_EZH1_CS_CKDL240001572-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_WT_EZH1_CS_CKDL240001572-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_WT_EZH1_CS_L14_2.fq.gz
cat input_raw_Novogene/NPC_WT_EZH2_CKDL240001573-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_WT_EZH2_CKDL240001573-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_WT_EZH2_L14_1.fq.gz
cat input_raw_Novogene/NPC_WT_EZH2_CKDL240001573-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_WT_EZH2_CKDL240001573-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_WT_EZH2_L14_2.fq.gz
cat input_raw_Novogene/NPC_WT_H3K27ac_CKDL240001571-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_WT_H3K27ac_CKDL240001571-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_WT_H3K27ac_L14_1.fq.gz
cat input_raw_Novogene/NPC_WT_H3K27ac_CKDL240001571-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_WT_H3K27ac_CKDL240001571-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_WT_H3K27ac_L14_2.fq.gz
cat input_raw_Novogene/NPC_WT_H3K4me3_CKDL240001569-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_WT_H3K4me3_CKDL240001569-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_WT_H3K4me3_L14_1.fq.gz
cat input_raw_Novogene/NPC_WT_H3K4me3_CKDL240001569-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_WT_H3K4me3_CKDL240001569-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_WT_H3K4me3_L14_2.fq.gz
cat input_raw_Novogene/NPC_WT_SUZ12_CKDL240001574-1A_H2NMCDSXC_L1_1.fq.gz input_raw_Novogene/NPC_WT_SUZ12_CKDL240001574-1A_HWK3JDSX7_L4_1.fq.gz > input/NPC_WT_SUZ12_L14_1.fq.gz
cat input_raw_Novogene/NPC_WT_SUZ12_CKDL240001574-1A_H2NMCDSXC_L1_2.fq.gz input_raw_Novogene/NPC_WT_SUZ12_CKDL240001574-1A_HWK3JDSX7_L4_2.fq.gz > input/NPC_WT_SUZ12_L14_2.fq.gz



