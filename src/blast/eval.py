import numpy as np
from sklearn.metrics import roc_auc_score


def iou(correct_range, matched_range):
    intersect = max(0, min(correct_range[1], matched_range[1]) - max(correct_range[0], matched_range[0]) + 1)
    union = max(correct_range[1], matched_range[1]) - min(correct_range[0], matched_range[0]) + 1
    iou = intersect / union
    return iou


def compute_performance(correct_range, outputs):

    label_pos = correct_range[0]
    for out in outputs:
        matched_range = [out['final_result']['ref_left_idx'], out['final_result']['ref_right_idx']]
        out['iou'] = iou(correct_range, matched_range)
        seed_matching_result = out['seed_matching_result']
        out['seed_matched'] = seed_matching_result['query_idx'] + label_pos == seed_matching_result['ref_idx']
        ungapped_extension_result = out['ungapped_extension_result']
        ungapped_range = [ungapped_extension_result['ref_left_idx'], ungapped_extension_result['ref_right_idx']]
        out['ungapped_extension_iou'] = iou(correct_range, ungapped_range)

    labels = [x['iou'] > 0 for x in outputs]
    scores = [x['final_result']['score'] for x in outputs]

    if len([x for x in outputs if x['iou'] > 0]) == 0:
        relevant_irelevant_auc = 0

    elif len([x for x in outputs if x['iou'] == 0]) == 0:
        relevant_irelevant_auc = 1
    else:
        relevant_irelevant_auc = roc_auc_score(labels, scores)

    top1_iou = outputs[0]['iou']
    top5_iou = np.max([x['iou'] for x in outputs[:5]])
    top10_iou = np.max([x['iou'] for x in outputs[:10]])

    seed_match_rate = len([x for x in outputs if x['seed_matched']]) / 90

    performance = {
        'relevant_irelevant_auc': relevant_irelevant_auc,
        'top1_iou': top1_iou,
        'top5_iou': top5_iou,
        'top10_iou': top10_iou,
        'seed_match_rate': seed_match_rate
    }
    return performance
