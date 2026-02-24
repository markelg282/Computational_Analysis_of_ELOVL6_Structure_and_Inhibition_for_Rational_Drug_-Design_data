import mdtraj as md
from pathlib import Path


def extract(structure_file, trajectory_file, num_snapshots, output_dir):
    
    traj = md.load(trajectory_file, top=structure_file)
    
    num_frames = len(traj)
    last_processed_frame = num_frames - (num_frames % num_snapshots)
    sample_freq = num_frames // num_snapshots
    sampled_frames = traj[0:last_processed_frame:sample_freq]
    
    for i, frame in enumerate(sampled_frames):
        frame_file = Path(output_dir, f'structure_{i+1}.gro')
        frame.save_gro(frame_file.as_posix())


if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-s', '--struct', type=str, required=True, help='Structure file for the pulling trajectory.')
    parser.add_argument('-t', '--traj', type=str, required=True, help='Trajectory file for the pulling trajectory.')
    parser.add_argument('-n', '--num_frames', type=int, required=True, help='Number of frames to extract from the pulling simulation.')
    parser.add_argument('-o', '--out_dir', type=str, required=True, help='Directory where output frame structures will be stored.')
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    
    extract(args.struct, args.traj, args.num_frames, out_dir.as_posix())
