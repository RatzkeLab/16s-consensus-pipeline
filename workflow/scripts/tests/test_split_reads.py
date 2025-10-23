#!/usr/bin/env python3
"""
Unit tests for split_reads.py
"""
import sys
from pathlib import Path
import tempfile

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from split_reads import read_cluster_assignments, parse_fastq


def test_read_cluster_assignments():
    """Test reading cluster assignments TSV with various header formats."""
    print("\n→ Testing read_cluster_assignments()...")
    
    # Create temporary TSV with ONT-style headers
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write("read_id\tcluster\n")
        f.write("25c35a0f-a899-4d31-a1f2-be93d5d70727 runid=abc ch=1\tA\n")
        f.write("3f4e5d6c-1234-5678-9abc-def012345678 runid=xyz ch=2\tB\n")
        f.write("12345678-aaaa-bbbb-cccc-dddddddddddd runid=xyz ch=3\tA\n")
        tsv_path = f.name
    
    try:
        assignments = read_cluster_assignments(tsv_path)
        
        # Check all reads are present with just UUID as key
        assert len(assignments) == 3, f"Expected 3 assignments, got {len(assignments)}"
        assert "25c35a0f-a899-4d31-a1f2-be93d5d70727" in assignments, "First UUID not found"
        assert "3f4e5d6c-1234-5678-9abc-def012345678" in assignments, "Second UUID not found"
        assert "12345678-aaaa-bbbb-cccc-dddddddddddd" in assignments, "Third UUID not found"
        
        # Check cluster assignments
        assert assignments["25c35a0f-a899-4d31-a1f2-be93d5d70727"] == "A"
        assert assignments["3f4e5d6c-1234-5678-9abc-def012345678"] == "B"
        assert assignments["12345678-aaaa-bbbb-cccc-dddddddddddd"] == "A"
        
        print("  ✓ Correctly extracts UUID from full headers")
        print("  ✓ Correctly assigns clusters")
    finally:
        Path(tsv_path).unlink()


def test_parse_fastq():
    """Test FASTQ parsing with ONT-style headers."""
    print("\n→ Testing parse_fastq()...")
    
    # Create temporary FASTQ with ONT-style headers
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write("@25c35a0f-a899-4d31-a1f2-be93d5d70727 runid=abc ch=1 start_time=2025\n")
        f.write("ACGTACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIII\n")
        f.write("@3f4e5d6c-1234-5678-9abc-def012345678 runid=xyz ch=2 start_time=2025\n")
        f.write("TGCATGCATGCATGCA\n")
        f.write("+\n")
        f.write("JJJJJJJJJJJJJJJJ\n")
        fastq_path = f.name
    
    try:
        reads = list(parse_fastq(fastq_path))
        
        assert len(reads) == 2, f"Expected 2 reads, got {len(reads)}"
        
        # Check first read
        read_id1, header1, seq1, qual1 = reads[0]
        assert read_id1 == "25c35a0f-a899-4d31-a1f2-be93d5d70727", f"Wrong read ID: {read_id1}"
        assert header1.startswith("@25c35a0f"), f"Wrong header: {header1}"
        assert seq1 == "ACGTACGTACGTACGT", f"Wrong sequence: {seq1}"
        assert qual1 == "IIIIIIIIIIIIIIII", f"Wrong quality: {qual1}"
        
        # Check second read
        read_id2, header2, seq2, qual2 = reads[1]
        assert read_id2 == "3f4e5d6c-1234-5678-9abc-def012345678", f"Wrong read ID: {read_id2}"
        assert seq2 == "TGCATGCATGCATGCA", f"Wrong sequence: {seq2}"
        
        print("  ✓ Correctly extracts UUID from FASTQ headers")
        print("  ✓ Preserves full headers, sequences, and quality scores")
    finally:
        Path(fastq_path).unlink()


def test_full_split_workflow():
    """Test the complete split workflow with realistic ONT data."""
    print("\n→ Testing full split workflow...")
    
    # Create cluster assignments with full ONT-style headers
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write("read_id\tcluster\n")
        f.write("25c35a0f-a899-4d31-a1f2-be93d5d70727 runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=440 start_time=2025-08-08T18:56:00.509097+02:00 flow_cell_id=FBC62275\tA\n")
        f.write("3f4e5d6c-1234-5678-9abc-def012345678 runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=100 start_time=2025-08-08T19:00:00.509097+02:00 flow_cell_id=FBC62275\tB\n")
        f.write("aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=200 start_time=2025-08-08T19:05:00.509097+02:00 flow_cell_id=FBC62275\tA\n")
        tsv_path = f.name
    
    # Create FASTQ with matching headers (as they appear in actual data)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write("@25c35a0f-a899-4d31-a1f2-be93d5d70727 runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=440 start_time=2025-08-08T18:56:00.509097+02:00 flow_cell_id=FBC62275 basecall_gpu=NVIDIA parent_read_id=25c35a0f-a899-4d31-a1f2-be93d5d70727\n")
        f.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
        f.write("@3f4e5d6c-1234-5678-9abc-def012345678 runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=100 start_time=2025-08-08T19:00:00.509097+02:00 flow_cell_id=FBC62275 basecall_gpu=NVIDIA parent_read_id=3f4e5d6c-1234-5678-9abc-def012345678\n")
        f.write("TGCATGCATGCATGCATGCATGCATGCATGCA\n")
        f.write("+\n")
        f.write("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ\n")
        f.write("@aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee runid=2ea20416c531223c9012fdc2b29405388434ef48 ch=200 start_time=2025-08-08T19:05:00.509097+02:00 flow_cell_id=FBC62275 basecall_gpu=NVIDIA parent_read_id=aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee\n")
        f.write("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        f.write("+\n")
        f.write("KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK\n")
        f.write("@unassigned-read-0000-0000-000000000000 runid=abc ch=999\n")
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n")
        f.write("+\n")
        f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        fastq_path = f.name
    
    try:
        # Read assignments
        assignments = read_cluster_assignments(tsv_path)
        print(f"  Loaded {len(assignments)} cluster assignments")
        
        # Create output directory
        with tempfile.TemporaryDirectory() as outdir:
            outdir_path = Path(outdir)
            
            # Open output files
            cluster_files = {}
            clusters_seen = set(assignments.values())
            
            for cluster in sorted(clusters_seen):
                output_path = outdir_path / f"test_{cluster}.fastq"
                cluster_files[cluster] = open(output_path, 'w')
            
            # Split reads
            reads_written = {cluster: 0 for cluster in clusters_seen}
            reads_unassigned = 0
            
            for read_id, header, seq, qual in parse_fastq(fastq_path):
                if read_id in assignments:
                    cluster = assignments[read_id]
                    cluster_files[cluster].write(f"{header}\n{seq}\n+\n{qual}\n")
                    reads_written[cluster] += 1
                else:
                    reads_unassigned += 1
            
            # Close files
            for f in cluster_files.values():
                f.close()
            
            # Verify results
            assert reads_written["A"] == 2, f"Expected 2 reads in cluster A, got {reads_written['A']}"
            assert reads_written["B"] == 1, f"Expected 1 read in cluster B, got {reads_written['B']}"
            assert reads_unassigned == 1, f"Expected 1 unassigned read, got {reads_unassigned}"
            
            print(f"  ✓ Cluster A: {reads_written['A']} reads")
            print(f"  ✓ Cluster B: {reads_written['B']} reads")
            print(f"  ✓ Unassigned: {reads_unassigned} reads")
            
            # Verify cluster A content
            cluster_a = outdir_path / "test_A.fastq"
            with open(cluster_a) as f:
                content_a = f.read()
                assert "25c35a0f-a899-4d31-a1f2-be93d5d70727" in content_a
                assert "ACGTACGTACGTACGTACGTACGTACGTACGT" in content_a
                assert "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee" in content_a
                assert "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" in content_a
                # Make sure cluster B read is NOT in cluster A
                assert "3f4e5d6c-1234-5678-9abc-def012345678" not in content_a
            
            # Verify cluster B content
            cluster_b = outdir_path / "test_B.fastq"
            with open(cluster_b) as f:
                content_b = f.read()
                assert "3f4e5d6c-1234-5678-9abc-def012345678" in content_b
                assert "TGCATGCATGCATGCATGCATGCATGCATGCA" in content_b
                # Make sure cluster A reads are NOT in cluster B
                assert "25c35a0f-a899-4d31-a1f2-be93d5d70727" not in content_b
                assert "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee" not in content_b
            
            print("  ✓ Cluster files contain correct reads")
            print("  ✓ No cross-contamination between clusters")
    finally:
        Path(tsv_path).unlink()
        Path(fastq_path).unlink()


if __name__ == "__main__":
    print("=" * 60)
    print("Testing split_reads.py")
    print("=" * 60)
    
    test_read_cluster_assignments()
    test_parse_fastq()
    test_full_split_workflow()
    
    print("\n" + "=" * 60)
    print("✓ All tests passed!")
    print("=" * 60)
