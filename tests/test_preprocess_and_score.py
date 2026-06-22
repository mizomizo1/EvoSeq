import tempfile
import unittest
from pathlib import Path

from evoseq.preprocess import preprocess_files, preprocess_folder
from evoseq.scoring import score_pairs_file
from evoseq.scoring.evo2_model import _import_evo2_class


class FakeScorer:
    def score_sequences(self, sequences, batch_size=8, label="sequences", progress=True):
        return [float(len(seq)) for seq in sequences]


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


class EvoSeqWorkflowTests(unittest.TestCase):
    def make_input_dir(self):
        tmp = tempfile.TemporaryDirectory()
        root = Path(tmp.name)
        input_dir = root / "test"

        write_text(
            input_dir / "evo2_reference.fasta",
            ">src|s1|GENE1|locus1|A>T\nACGTACGT\n"
            ">src|s2|GENE2|locus2|G>C\nGGGGCCCC\n",
        )
        write_text(
            input_dir / "evo2_mutant.fasta",
            ">src|s1|GENE1|locus1|A>T\nACGTTCGT\n"
            ">src|s2|GENE2|locus2|G>C\nGGGGACCC\n",
        )
        write_text(
            input_dir / "manifest.tsv",
            "record_id\tsample\tgene\thgvs\tannotation\n"
            "src|s1|GENE1|locus1|A>T\tS1\tGENE1\tNM_1:c.1A>T\tmissense\n"
            "src|s2|GENE2|locus2|G>C\tS2\tGENE2\tNM_2:c.2G>C\tsplicing\n",
        )
        return tmp, input_dir

    def test_preprocess_files_writes_dedicated_output_dir(self):
        tmp, input_dir = self.make_input_dir()
        self.addCleanup(tmp.cleanup)

        table, paths = preprocess_files(
            reference_fasta_path=input_dir / "evo2_reference.fasta",
            mutant_fasta_path=input_dir / "evo2_mutant.fasta",
            manifest_path="auto",
            progress=False,
        )

        self.assertEqual(len(table), 2)
        self.assertEqual(
            paths["pairs"].parent.resolve(),
            (input_dir / "evoseq_preprocess_output").resolve(),
        )
        self.assertTrue(paths["pairs"].exists())
        self.assertTrue(paths["preprocessing_report"].exists())
        self.assertIn("annotation", table.columns)

    def test_preprocess_folder_discovers_files(self):
        tmp, input_dir = self.make_input_dir()
        self.addCleanup(tmp.cleanup)

        table, paths = preprocess_folder(input_dir, progress=False)

        self.assertEqual(len(table), 2)
        self.assertTrue(paths["reference"].exists())
        self.assertTrue(paths["mutant"].exists())

    def test_score_pairs_file_with_fake_scorer_does_not_need_evo2_or_torch(self):
        tmp, input_dir = self.make_input_dir()
        self.addCleanup(tmp.cleanup)

        _, preprocess_paths = preprocess_files(
            reference_fasta_path=input_dir / "evo2_reference.fasta",
            mutant_fasta_path=input_dir / "evo2_mutant.fasta",
            manifest_path="auto",
            progress=False,
        )

        result, score_paths = score_pairs_file(
            preprocess_paths["pairs"],
            scorer=FakeScorer(),
            device="cpu",
            progress=False,
        )

        self.assertEqual(len(result), 2)
        self.assertTrue(score_paths["unique_scores"].exists())
        self.assertTrue(score_paths["scoring_report"].exists())
        self.assertEqual(
            score_paths["unique_scores"].parent.resolve(),
            (preprocess_paths["pairs"].parent / "evoseq_scoring_output").resolve(),
        )
        self.assertIn("evo2_delta_score", result.columns)

    def test_missing_evo2_error_mentions_runtime_dependencies(self):
        try:
            _import_evo2_class()
        except RuntimeError as exc:
            message = str(exc)
            self.assertIn("Evo2", message)
            self.assertIn("flash-attn", message)
            self.assertIn("torch", message)


if __name__ == "__main__":
    unittest.main()
