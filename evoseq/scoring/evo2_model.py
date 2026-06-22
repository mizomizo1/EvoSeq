import torch


def load_evo2_model(model_name="evo2_7b", device="cuda:0"):
    from evo2 import Evo2

    model = Evo2(model_name)
    return model, device


def tokenize_sequence(model, seq, device):
    token_ids = model.tokenizer.tokenize(seq)
    return torch.tensor(token_ids, dtype=torch.int, device=device).unsqueeze(0)