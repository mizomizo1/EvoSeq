import torch
import torch.nn.functional as F

from .evo2_model import tokenize_sequence


def per_token_logprob(model, seq, device="cuda:0"):
    """
    Return per-position log-probability for a DNA sequence.

    Position 0 has no preceding context, so its log probability is None.
    """

    if not seq:
        return []

    ids = tokenize_sequence(model, seq, device)

    with torch.inference_mode():
        outputs, _ = model(ids)
        logits = outputs[0] if isinstance(outputs, (tuple, list)) else outputs

    if logits.ndim != 3:
        raise RuntimeError(
            f"Expected logits [batch, length, vocab], got {tuple(logits.shape)}"
        )

    if tuple(logits.shape[:2]) != tuple(ids.shape):
        raise RuntimeError(
            f"Logits/input shape mismatch: logits={tuple(logits.shape)}, ids={tuple(ids.shape)}"
        )

    if ids.shape[1] != len(seq):
        raise RuntimeError(
            f"Token/base length mismatch: tokens={ids.shape[1]}, bases={len(seq)}"
        )

    shifted_logprobs = F.log_softmax(logits[:, :-1, :].float(), dim=-1)
    target_ids = ids[:, 1:].long().unsqueeze(-1)

    per_position = shifted_logprobs.gather(
        dim=-1,
        index=target_ids,
    ).squeeze(-1)[0]

    per_position = per_position.cpu().tolist()

    return [(seq[0], None)] + [
        (seq[i], per_position[i - 1])
        for i in range(1, len(seq))
    ]


def extract_center_window(per_token_results, center=4096, half_window=320):
    rows = []

    lo = center - half_window
    hi = center + half_window

    for i in range(max(1, lo), min(len(per_token_results), hi)):
        base, logprob = per_token_results[i]

        if logprob is None:
            continue

        rows.append(
            {
                "pos_index_0based": i,
                "rel_to_center": i - center,
                "base": base,
                "logprob": logprob,
            }
        )

    return rows