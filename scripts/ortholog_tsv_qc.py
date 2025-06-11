import pandas as pd
import csv
import sys
from typing import Dict, List, Tuple, Set

def qc_tsv_file(input_tsv: str, expected_columns: List[str] = None) -> Dict:
    """
    Quality control check for TSV file to identify malformed lines and potential data loss.
    
    Args:
        input_tsv: Path to input TSV file
        expected_columns: List of expected column names
        
    Returns:
        Dictionary with QC results and statistics
    """
    
    if expected_columns is None:
        expected_columns = [
            'query', 'target', 'query_id', 'scenario', 'query_contig', 'query_start', 
            'query_end', 'target_contig', 'target_start', 'target_end', 'query_gene', 
            'target_gene', 'query_family', 'target_family', 'query_splice_site', 
            'target_splice_site', 'query_introner_id', 'target_introner_id', 
            'match_score', 'orientation', 'left_reverse', 'right_reverse'
        ]
    
    qc_results = {
        'file_path': input_tsv,
        'total_lines': 0,
        'header_line': None,
        'data_lines': 0,
        'malformed_lines': [],
        'column_count_issues': [],
        'missing_expected_columns': [],
        'extra_columns': [],
        'empty_lines': 0,
        'pandas_read_successful': True,
        'pandas_row_count': 0,
        'row_count_mismatch': False,
        'data_type_issues': [],
        'coordinate_issues': []
    }
    
    print(f"Starting QC analysis of {input_tsv}")
    print("="*60)
    
    # Step 1: Raw line-by-line analysis
    print("Step 1: Analyzing raw file structure...")
    
    try:
        with open(input_tsv, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            qc_results['total_lines'] = len(lines)
            
            if not lines:
                print("ERROR: File is empty!")
                return qc_results
            
            # Check header
            header_line = lines[0].strip()
            qc_results['header_line'] = header_line
            header_cols = header_line.split('\t')
            
            print(f"Header found with {len(header_cols)} columns")
            print(f"Expected {len(expected_columns)} columns")
            
            # Check for missing/extra columns
            header_set = set(header_cols)
            expected_set = set(expected_columns)
            
            qc_results['missing_expected_columns'] = list(expected_set - header_set)
            qc_results['extra_columns'] = list(header_set - expected_set)
            
            if qc_results['missing_expected_columns']:
                print(f"WARNING: Missing expected columns: {qc_results['missing_expected_columns']}")
            
            if qc_results['extra_columns']:
                print(f"INFO: Extra columns found: {qc_results['extra_columns']}")
            
            # Analyze each data line
            expected_col_count = len(header_cols)
            column_count_distribution = {}  # Track distribution of column counts
            malformed_patterns = {}  # Track patterns in malformed lines
            
            for line_num, line in enumerate(lines[1:], start=2):  # Start from line 2 (after header)
                line_stripped = line.strip()
                
                if not line_stripped:
                    qc_results['empty_lines'] += 1
                    continue
                
                qc_results['data_lines'] += 1
                cols = line_stripped.split('\t')
                actual_col_count = len(cols)
                
                # Track column count distribution
                if actual_col_count not in column_count_distribution:
                    column_count_distribution[actual_col_count] = 0
                column_count_distribution[actual_col_count] += 1
                
                # Check column count
                if actual_col_count != expected_col_count:
                    # Track patterns in malformed lines
                    pattern_key = f"{actual_col_count}_cols"
                    if pattern_key not in malformed_patterns:
                        malformed_patterns[pattern_key] = {
                            'count': 0,
                            'examples': [],
                            'first_occurrence': line_num
                        }
                    
                    malformed_patterns[pattern_key]['count'] += 1
                    
                    # Store first few examples
                    if len(malformed_patterns[pattern_key]['examples']) < 3:
                        malformed_patterns[pattern_key]['examples'].append({
                            'line_number': line_num,
                            'line_preview': line_stripped[:150] + '...' if len(line_stripped) > 150 else line_stripped
                        })
                    
                    qc_results['malformed_lines'].append(line_num)
            
            # Add analysis results
            qc_results['column_count_distribution'] = column_count_distribution
            qc_results['malformed_patterns'] = malformed_patterns
    
    except Exception as e:
        print(f"ERROR reading file: {e}")
        qc_results['file_read_error'] = str(e)
        return qc_results
    
    print(f"Raw analysis complete:")
    print(f"  - Total lines: {qc_results['total_lines']}")
    print(f"  - Data lines: {qc_results['data_lines']}")
    print(f"  - Empty lines: {qc_results['empty_lines']}")
    print(f"  - Malformed lines: {len(qc_results['malformed_lines'])}")
    
    # Print column count distribution
    if 'column_count_distribution' in qc_results:
        print(f"\nColumn Count Distribution:")
        for col_count, frequency in sorted(qc_results['column_count_distribution'].items()):
            status = "✅" if col_count == len(header_cols) else "❌"
            print(f"  {status} {col_count} columns: {frequency:,} lines ({frequency/qc_results['data_lines']*100:.1f}%)")
    
    # Print malformed patterns
    if 'malformed_patterns' in qc_results and qc_results['malformed_patterns']:
        print(f"\nMalformed Line Patterns:")
        for pattern, details in sorted(qc_results['malformed_patterns'].items(), 
                                     key=lambda x: x[1]['count'], reverse=True):
            col_count = int(pattern.split('_')[0])
            diff = col_count - len(header_cols)
            diff_str = f"({diff:+d})" if diff != 0 else ""
            print(f"  ❌ {col_count} columns {diff_str}: {details['count']:,} lines")
            print(f"     First seen at line {details['first_occurrence']}")
            if details['examples']:
                print(f"     Example: {details['examples'][0]['line_preview']}")
            print()
    
    # Step 2: Try pandas read and compare
    print("\nStep 2: Testing pandas read...")
    
    try:
        # Read with pandas (your current method)
        df = pd.read_csv(input_tsv, sep='\t')
        qc_results['pandas_row_count'] = len(df)
        qc_results['pandas_columns'] = list(df.columns)
        
        print(f"Pandas successfully read {len(df)} rows")
        
        # Check for row count mismatch
        if qc_results['pandas_row_count'] != qc_results['data_lines']:
            qc_results['row_count_mismatch'] = True
            dropped_rows = qc_results['data_lines'] - qc_results['pandas_row_count']
            print(f"WARNING: Row count mismatch! {dropped_rows} rows may have been dropped")
            print(f"  Raw file data lines: {qc_results['data_lines']}")
            print(f"  Pandas DataFrame rows: {qc_results['pandas_row_count']}")
        
        # Step 3: Data type and coordinate validation
        print("\nStep 3: Validating data types and coordinates...")
        
        # Check numeric columns that should be integers
        numeric_cols = ['query_start', 'query_end', 'target_start', 'target_end', 'scenario']
        for col in numeric_cols:
            if col in df.columns:
                # Check for non-numeric values (excluding NaN)
                non_numeric = df[col].apply(lambda x: pd.notna(x) and not str(x).replace('.', '').replace('-', '').isdigit())
                if non_numeric.any():
                    problematic_values = df.loc[non_numeric, col].unique()
                    qc_results['data_type_issues'].append({
                        'column': col,
                        'problematic_values': list(problematic_values),
                        'count': non_numeric.sum()
                    })
        
        # Check coordinate consistency
        coord_issues = []
        for idx, row in df.iterrows():
            # Check if start <= end for both query and target
            try:
                if (pd.notna(row.get('query_start')) and pd.notna(row.get('query_end')) and 
                    int(float(row['query_start'])) > int(float(row['query_end']))):
                    coord_issues.append({
                        'row': idx + 2,  # +2 because pandas is 0-indexed and we skip header
                        'issue': 'query_start > query_end',
                        'query_start': row['query_start'],
                        'query_end': row['query_end']
                    })
                
                if (pd.notna(row.get('target_start')) and pd.notna(row.get('target_end')) and 
                    int(float(row['target_start'])) > int(float(row['target_end']))):
                    coord_issues.append({
                        'row': idx + 2,
                        'issue': 'target_start > target_end', 
                        'target_start': row['target_start'],
                        'target_end': row['target_end']
                    })
            except (ValueError, TypeError):
                continue  # Skip rows with invalid coordinate data
        
        qc_results['coordinate_issues'] = coord_issues
        
    except Exception as e:
        print(f"ERROR: Pandas failed to read file: {e}")
        qc_results['pandas_read_successful'] = False
        qc_results['pandas_error'] = str(e)
    
    return qc_results

def print_qc_summary(qc_results: Dict):
    """Print a comprehensive QC summary"""
    
    print("\n" + "="*60)
    print("QC SUMMARY")
    print("="*60)
    
    # Overall status
    issues_found = (
        len(qc_results['malformed_lines']) > 0 or
        qc_results['row_count_mismatch'] or
        len(qc_results['missing_expected_columns']) > 0 or
        len(qc_results['data_type_issues']) > 0 or
        len(qc_results['coordinate_issues']) > 0
    )
    
    status = "❌ ISSUES FOUND" if issues_found else "✅ PASSED"
    print(f"Overall Status: {status}")
    print()
    
    # Detailed results
    if qc_results['malformed_lines']:
        print(f"❌ Malformed Lines: {len(qc_results['malformed_lines']):,}")
        
        # Show patterns instead of individual lines
        if 'malformed_patterns' in qc_results:
            print(f"   Breakdown by column count:")
            for pattern, details in sorted(qc_results['malformed_patterns'].items(), 
                                         key=lambda x: x[1]['count'], reverse=True):
                col_count = int(pattern.split('_')[0])
                expected = len(qc_results.get('pandas_columns', qc_results['header_line'].split('\t')))
                diff = col_count - expected
                diff_str = f" ({diff:+d} columns)" if diff != 0 else ""
                print(f"     • {col_count} columns{diff_str}: {details['count']:,} lines")
        print()
    
    if qc_results['row_count_mismatch']:
        dropped = qc_results['data_lines'] - qc_results['pandas_row_count']
        print(f"❌ Row Count Mismatch: {dropped} rows dropped by pandas")
        print(f"   Raw file: {qc_results['data_lines']} data lines")
        print(f"   Pandas read: {qc_results['pandas_row_count']} rows")
        print()
    
    if qc_results['missing_expected_columns']:
        print(f"❌ Missing Expected Columns: {qc_results['missing_expected_columns']}")
        print()
    
    if qc_results['data_type_issues']:
        print(f"❌ Data Type Issues:")
        for issue in qc_results['data_type_issues']:
            print(f"   Column '{issue['column']}': {issue['count']} non-numeric values")
            print(f"   Examples: {issue['problematic_values'][:5]}")
        print()
    
    if qc_results['coordinate_issues']:
        print(f"❌ Coordinate Issues: {len(qc_results['coordinate_issues'])}")
        for issue in qc_results['coordinate_issues'][:5]:  # Show first 5
            print(f"   Row {issue['row']}: {issue['issue']}")
        if len(qc_results['coordinate_issues']) > 5:
            print(f"   ... and {len(qc_results['coordinate_issues']) - 5} more")
        print()
    
    # Positive results
    if not issues_found:
        print("✅ File structure appears correct")
        print("✅ All expected columns present")
        print("✅ No malformed lines detected")
        print("✅ Row counts match between raw file and pandas")

def main():
    """Main function to run QC check"""
    
    if len(sys.argv) != 2:
        print("Usage: python qc_script.py <input_tsv_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Run QC
    results = qc_tsv_file(input_file)
    
    # Print summary
    print_qc_summary(results)
    
    # Return exit code based on results
    if (len(results['malformed_lines']) > 0 or 
        results['row_count_mismatch'] or 
        len(results['missing_expected_columns']) > 0):
        sys.exit(1)  # Exit with error code
    else:
        sys.exit(0)  # Success

if __name__ == "__main__":
    main()