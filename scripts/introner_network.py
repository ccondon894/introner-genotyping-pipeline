import pandas as pd
import networkx as nx
from collections import defaultdict
import os
import sys

def build_ortholog_groups(input_tsv):
    """
    Build ortholog groups from query-target pairs.
    
    Args:
        input_tsv: Path to input TSV from previous script
        
    Returns:
        DataFrame with ortholog groups
    """
    print(f"Reading input data from {input_tsv}...")
    
    # Read the input TSV
    df = pd.read_csv(input_tsv, sep='\t')
    print(f"Read {len(df)} query-target pairs")
    
    # DEBUG - Count distinct introners per sample in input file
    input_introners_by_sample = {}
    for sample in df['query'].unique():
        # Get all introners (unique coordinates) for this sample
        sample_rows = df[df['query'] == sample]
        unique_coords = set()
        for _, row in sample_rows.iterrows():
            if pd.notna(row['query_contig']) and pd.notna(row['query_start']) and pd.notna(row['query_end']):
                key = f"{row['query_contig']}:{row['query_start']}-{row['query_end']}"
                unique_coords.add(key)
        
        input_introners_by_sample[sample] = unique_coords
        print(f"Sample {sample} has {len(unique_coords)} distinct introners in input")
    
    # Filter out rows with scenario 3 (missing data)
    valid_df = df[df['scenario'] != 3].copy()
    print(f"After filtering scenario 3: {len(valid_df)} pairs")
    
    # DEBUG - Count how many introners were filtered due to scenario 3
    scenario3_filtered = {}
    for sample in df['query'].unique():
        # Get all introners for this sample in original df
        all_introners = input_introners_by_sample[sample]
        
        # Get remaining introners after filtering
        remaining_introners = set()
        for _, row in valid_df[valid_df['query'] == sample].iterrows():
            if pd.notna(row['query_contig']) and pd.notna(row['query_start']) and pd.notna(row['query_end']):
                key = f"{row['query_contig']}:{row['query_start']}-{row['query_end']}"
                remaining_introners.add(key)
        
        filtered_introners = all_introners - remaining_introners
        scenario3_filtered[sample] = filtered_introners
        print(f"Sample {sample} had {len(filtered_introners)} introners filtered due to scenario 3")
    
    # Create sample-specific identifiers for nodes
    # Explicitly convert coordinates to integers for the node IDs to prevent float formatting
    valid_df['query_node'] = valid_df['query'] + '|' + valid_df['query_contig'].astype(str) + ':' + \
                            valid_df['query_start'].apply(lambda x: str(int(x)) if pd.notnull(x) else '').astype(str) + '-' + \
                            valid_df['query_end'].apply(lambda x: str(int(x)) if pd.notnull(x) else '').astype(str)
    valid_df['target_node'] = valid_df['target'] + '|' + valid_df['target_contig'].astype(str) + ':' + \
                             valid_df['target_start'].apply(lambda x: str(int(x)) if pd.notnull(x) else '').astype(str) + '-' + \
                             valid_df['target_end'].apply(lambda x: str(int(x)) if pd.notnull(x) else '').astype(str)
    
    # Build a graph of relationships
    G = nx.Graph()
    
    # Add edges for each valid query-target pair
    for _, row in valid_df.iterrows():
        # Add nodes with metadata
        G.add_node(row['query_node'], 
                   sample=row['query'],
                   contig=row['query_contig'],
                   start=row['query_start'],
                   end=row['query_end'],
                   gene=row['query_gene'],
                   family=row['query_family'],
                   splice_site=row['query_splice_site'],
                   introner_id=row['query_introner_id'],
                   presence=1)  # Query always has the introner
                   
        G.add_node(row['target_node'], 
                   sample=row['target'],
                   contig=row['target_contig'],
                   start=row['target_start'],
                   end=row['target_end'],
                   gene=row['target_gene'],
                   family=row['target_family'],
                   splice_site=row['target_splice_site'],
                   introner_id=row['target_introner_id'],
                   presence=1 if row['scenario'] == 1 else 2)  # Target presence depends on scenario
                   
        # Add edge with match score
        G.add_edge(row['query_node'], row['target_node'], 
                   score=row.get('match_score', 0),
                   scenario=row['scenario'])
    
    print(f"Built graph with {len(G.nodes)} nodes and {len(G.edges)} edges")
    
    # DEBUG - Track which introners made it into the graph
    graph_introners_by_sample = defaultdict(set)
    for node in G.nodes:
        attrs = G.nodes[node]
        sample = attrs['sample']
        contig = attrs['contig']
        start = attrs['start']
        end = attrs['end']
        
        if pd.notna(contig) and pd.notna(start) and pd.notna(end):
            key = f"{contig}:{start}-{end}"
            graph_introners_by_sample[sample].add(key)
    
    # DEBUG - Check specific problematic introners
    problematic_introners = [
        "RCC3052|RCC3052-deplete#0#7#1:163243-163619",
        "RCC1749|RCC1749#0#intronerless_contig_25:163216-163592"
    ]
    
    print("\n*** DEBUG - Checking problematic introners ***")
    for introner_node in problematic_introners:
        if introner_node in G.nodes:
            print(f"\nFound problematic introner in graph: {introner_node}")
            print(f"Node attributes: {G.nodes[introner_node]}")
            neighbors = list(G.neighbors(introner_node))
            print(f"Number of neighbors: {len(neighbors)}")
            print("Neighbors:")
            for neighbor in neighbors:
                neighbor_attrs = G.nodes[neighbor]
                edge_attrs = G.edges[introner_node, neighbor]
                print(f"  - {neighbor}")
                print(f"    Sample: {neighbor_attrs['sample']}")
                print(f"    Contig: {neighbor_attrs['contig']}")
                print(f"    Position: {neighbor_attrs['start']}-{neighbor_attrs['end']}")
                print(f"    Edge score: {edge_attrs.get('score', 0)}")
                print(f"    Edge scenario: {edge_attrs.get('scenario', 'unknown')}")
    
    # DEBUG - Compare input introners to those in the graph
    print("\nComparison of introners in input vs. graph:")
    for sample in input_introners_by_sample:
        input_count = len(input_introners_by_sample[sample])
        graph_count = len(graph_introners_by_sample[sample])
        
        # Calculate loss at this stage
        lost_count = input_count - len(scenario3_filtered[sample]) - graph_count
        if lost_count > 0:
            print(f"Sample {sample}:")
            print(f"  Input: {input_count}")
            print(f"  Lost to scenario 3: {len(scenario3_filtered[sample])}")
            print(f"  In graph: {graph_count}")
            print(f"  Additional loss: {lost_count}")
            
            # Check if there are duplicates in the input that might explain the loss
            query_coords = [f"{row['query_contig']}:{row['query_start']}-{row['query_end']}" 
                            for _, row in valid_df[valid_df['query'] == sample].iterrows()
                            if pd.notna(row['query_contig']) and pd.notna(row['query_start']) and pd.notna(row['query_end'])]
            
            duplicate_count = len(query_coords) - len(set(query_coords))
            if duplicate_count > 0:
                print(f"  Found {duplicate_count} duplicate coordinates in input that might explain some loss")
    
    # Find connected components (potential ortholog groups)
    components = list(nx.connected_components(G))
    print(f"\nFound {len(components)} initial connected components")
    
    # DEBUG - Count component sizes for insight
    component_sizes = [len(comp) for comp in components]
    print(f"Component size stats - Min: {min(component_sizes)}, Max: {max(component_sizes)}, Avg: {sum(component_sizes)/len(component_sizes):.1f}")
    
    # Find component containing the problematic introners
    problematic_introners = [
        "RCC3052|RCC3052-deplete#0#7#1:163243-163619",
        "RCC1749|RCC1749#0#intronerless_contig_25:163216-163592"
    ]
    
    problematic_component = None
    for i, component in enumerate(components):
        for introner_node in problematic_introners:
            if introner_node in component:
                problematic_component = (i, component)
                break
        if problematic_component:
            break
    
    if problematic_component:
        comp_idx, component = problematic_component
        print(f"\nFound problematic introners in component {comp_idx} with {len(component)} nodes")
        
        # Count samples in this component
        samples_in_component = defaultdict(int)
        for node in component:
            sample = G.nodes[node]['sample']
            samples_in_component[sample] += 1
        
        print("Samples in component:")
        for sample, count in samples_in_component.items():
            print(f"  - {sample}: {count} nodes")
            
        # Check for any nodes from same introner_id but different samples
        introner_ids_to_nodes = defaultdict(list)
        for node in component:
            introner_id = G.nodes[node].get('introner_id', '')
            if introner_id:
                introner_ids_to_nodes[introner_id].append(node)
        
        print("\nIntroner IDs shared across different samples:")
        for introner_id, nodes in introner_ids_to_nodes.items():
            if len(nodes) > 1:
                samples = {G.nodes[node]['sample'] for node in nodes}
                if len(samples) > 1:
                    print(f"  - Introner ID {introner_id} found in samples: {', '.join(samples)}")
                    for node in nodes:
                        attrs = G.nodes[node]
                        print(f"    * {node} ({attrs['sample']}) at {attrs['contig']}:{attrs['start']}-{attrs['end']}")
    else:
        print("\nProblematic introners not found in any component")
    
    # Resolve conflicts where a sample appears multiple times in a component
    resolved_groups = []
    removed_nodes_count = 0
    removed_nodes_by_sample = defaultdict(int)
    
    for i, component in enumerate(components):
        # Create a map of sample to nodes
        sample_to_nodes = defaultdict(list)
        for node in component:
            sample = G.nodes[node]['sample']
            sample_to_nodes[sample].append(node)
        
        # Check if we have any conflicts (a sample appears multiple times)
        has_conflicts = any(len(nodes) > 1 for nodes in sample_to_nodes.values())
        
        if not has_conflicts:
            # No conflicts, add the component as is
            resolved_groups.append(list(component))
        else:
            # DEBUG - Count conflicts
            conflict_samples = [sample for sample, nodes in sample_to_nodes.items() if len(nodes) > 1]
            total_conflicts = sum(len(nodes) - 1 for sample, nodes in sample_to_nodes.items() if len(nodes) > 1)
            print(f"Component {i} has conflicts: {len(conflict_samples)} samples with {total_conflicts} total conflicts")
            
            # Resolve conflicts by creating a subgraph and finding the maximum weight clique
            # This is a simplification - in a real implementation you'd need a more sophisticated algorithm
            subgraph = G.subgraph(component).copy()
            component_removed_nodes = []
            
            # Keep removing lowest-weighted conflicting nodes until no conflicts remain
            while has_conflicts:
                # Find the sample with conflicts
                conflict_sample = next(sample for sample, nodes in sample_to_nodes.items() if len(nodes) > 1)
                conflict_nodes = sample_to_nodes[conflict_sample]
                
                # Find the node with the lowest total edge weight
                node_weights = {}
                for node in conflict_nodes:
                    node_weights[node] = sum(subgraph[node][neighbor].get('score', 0) for neighbor in subgraph[node])
                
                # Remove the lowest-weighted node
                lowest_node = min(node_weights, key=node_weights.get)
                component_removed_nodes.append((lowest_node, conflict_sample))
                subgraph.remove_node(lowest_node)
                removed_nodes_count += 1
                removed_nodes_by_sample[conflict_sample] += 1
                
                # Update sample_to_nodes
                sample_to_nodes = defaultdict(list)
                for node in subgraph.nodes:
                    sample = G.nodes[node]['sample']
                    sample_to_nodes[sample].append(node)
                
                # Check if we still have conflicts
                has_conflicts = any(len(nodes) > 1 for nodes in sample_to_nodes.values())
            
            # Add the resolved component
            resolved_groups.append(list(subgraph.nodes))
    
    print(f"\nAfter resolving conflicts: {len(resolved_groups)} ortholog groups")
    print(f"Removed {removed_nodes_count} nodes during conflict resolution")
    print("Nodes removed by sample during conflict resolution:")
    for sample, count in sorted(removed_nodes_by_sample.items()):
        print(f"  {sample}: {count}")
    
    # Create a mapping of nodes to ortholog groups
    node_to_group = {}
    for i, group in enumerate(resolved_groups):
        for node in group:
            node_to_group[node] = f"ortholog_id_{i+1:04d}"
    
    # Create the output dataframe
    rows = []
    for node, ortholog_id in node_to_group.items():
        # Get node attributes
        attrs = G.nodes[node]
        rows.append({
            'ortholog_id': ortholog_id,
            'sample': attrs['sample'],
            'start': attrs['start'],
            'end': attrs['end'],
            'contig': attrs['contig'],
            'sequence_id': node.split('|')[1] if '|' in node else '',
            'family': attrs.get('family', ''),
            'gene': attrs.get('gene', ''),
            'splice_site': attrs.get('splice_site', ''),
            'presence': attrs['presence']
        })
    
    # DEBUG - Track which introners made it to the final output
    final_introners_by_sample = defaultdict(set)
    for row in rows:
        sample = row['sample']
        contig = row['contig']
        start = row['start']
        end = row['end']
        
        if row['presence'] == 1 and pd.notna(contig) and pd.notna(start) and pd.notna(end):
            key = f"{contig}:{start}-{end}"
            final_introners_by_sample[sample].add(key)
    
    # DEBUG - Compare graph introners to final output
    print("\nComparison of introners in graph vs. final output:")
    for sample in input_introners_by_sample:
        if sample in graph_introners_by_sample and sample in final_introners_by_sample:
            graph_count = len(graph_introners_by_sample[sample])
            final_count = len(final_introners_by_sample[sample])
            lost_in_resolution = graph_count - final_count
            
            print(f"Sample {sample}:")
            print(f"  In graph: {graph_count}")
            print(f"  In final output: {final_count}")
            print(f"  Lost during conflict resolution: {lost_in_resolution}")
    
    # Add missing samples with scenario 3
    all_samples = set(df['query'].unique()) | set(df['target'].unique())
    groups_by_sample = defaultdict(list)
    
    for row in rows:
        groups_by_sample[row['sample']].append(row['ortholog_id'])
    
    for sample in all_samples:
        assigned_groups = groups_by_sample.get(sample, [])
        if not assigned_groups:
            # This sample isn't in any group, add it with scenario 3
            rows.append({
                'ortholog_id': f"ortholog_id_{len(resolved_groups)+1:04d}",
                'sample': sample,
                'start': '',
                'end': '',
                'contig': '',
                'sequence_id': '',
                'family': '',
                'gene': '',
                'splice_site': '',
                'presence': 3
            })
    
    # Convert to DataFrame
    result_df = pd.DataFrame(rows)
    
    # Check that each sample appears only once per ortholog group
    validation = result_df.groupby(['ortholog_id', 'sample']).size().reset_index(name='count')
    duplicates = validation[validation['count'] > 1]
    if not duplicates.empty:
        print(f"\nWARNING: Found {len(duplicates)} cases where a sample appears multiple times in an ortholog group")
    
    # DEBUG - Create a summary of introner loss at each stage
    loss_summary = []
    for sample in input_introners_by_sample:
        input_count = len(input_introners_by_sample[sample])
        scenario3_count = len(scenario3_filtered.get(sample, set()))
        graph_count = len(graph_introners_by_sample.get(sample, set()))
        final_count = len(final_introners_by_sample.get(sample, set()))
        removed_count = removed_nodes_by_sample.get(sample, 0)
        
        # Calculate total loss and percentage
        total_loss = input_count - final_count
        loss_percent = (total_loss / input_count) * 100 if input_count > 0 else 0
        
        loss_summary.append({
            'sample': sample,
            'input_introners': input_count,
            'lost_to_scenario3': scenario3_count,
            'added_to_graph': graph_count,
            'removed_in_conflict_resolution': removed_count,
            'final_output': final_count,
            'total_loss': total_loss,
            'loss_percent': loss_percent
        })
    
    
    return result_df

def fill_missing_samples(df, all_samples):
    """
    Fill in missing samples for each ortholog group with presence=3.
    
    Args:
        df: DataFrame with ortholog groups
        all_samples: List of all sample names
        
    Returns:
        DataFrame with missing samples added
    """
    # Get all ortholog IDs
    ortholog_ids = df['ortholog_id'].unique()
    
    # Create a set of (ortholog_id, sample) pairs that already exist
    existing_pairs = set(zip(df['ortholog_id'], df['sample']))
    
    # Create rows for missing samples
    missing_rows = []
    for oid in ortholog_ids:
        for sample in all_samples:
            if (oid, sample) not in existing_pairs:
                missing_rows.append({
                    'ortholog_id': oid,
                    'sample': sample,
                    'start': '',
                    'end': '',
                    'contig': '',
                    'sequence_id': '',
                    'family': '',
                    'gene': '',
                    'splice_site': '',
                    'presence': 3
                })
    
    # Append missing rows to the dataframe
    if missing_rows:
        missing_df = pd.DataFrame(missing_rows)
        df = pd.concat([df, missing_df], ignore_index=True)
    
    return df

def main():
    input_tsv = sys.argv[1]
    output_tsv = sys.argv[2]
    
    # Build ortholog groups
    result_df = build_ortholog_groups(input_tsv)
    
    # Get all unique samples
    all_samples = list(result_df['sample'].unique())
    print(f"Found {len(all_samples)} unique samples")
    
    # Fill in missing samples
    result_df = fill_missing_samples(result_df, all_samples)
    
    # Sort by ortholog_id and sample
    result_df = result_df.sort_values(['ortholog_id', 'sample'])
    
    # Convert numeric columns to the correct types
    # For start and end, convert to integers but preserve empty strings
    for col in ['start', 'end', 'family']:
        # First, check if the column exists
        if col in result_df.columns:
            # Convert to integers where possible, but preserve empty strings
            result_df[col] = result_df[col].apply(
                lambda x: int(x) if pd.notnull(x) and str(x).strip() and not pd.isna(x) else x
            )
    
    # Write to TSV with appropriate types
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Wrote {len(result_df)} rows to {output_tsv}")
    
    # Generate summary
    ortho_count = len(result_df['ortholog_id'].unique())
    presence_counts = result_df['presence'].value_counts()
    print(f"Generated {ortho_count} ortholog groups")
    print(f"Presence counts: {presence_counts.to_dict()}")
    
    # Generate final copy counts by sample
    final_counts = result_df[result_df['presence'] == 1].groupby('sample').size().reset_index(name='count')
    print("\nFinal present introner counts by sample:")
    for _, row in final_counts.iterrows():
        print(f"  {row['sample']}: {row['count']}")

if __name__ == "__main__":
    main()