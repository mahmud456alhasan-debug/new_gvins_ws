#!/usr/bin/env python3
"""Play a bag with paths that contain spaces (roslaunch splits plain args on whitespace)."""
import os
import subprocess
import sys


def main():
    bag = os.environ.get("REPLAY_BAG")
    if not bag:
        print("replay_bag.py: REPLAY_BAG environment variable not set", file=sys.stderr)
        return 1
    rate = os.environ.get("REPLAY_RATE", "1.0")
    start = os.environ.get("REPLAY_START", "0")
    distro = os.environ.get("ROS_DISTRO", "noetic")
    play = "/opt/ros/{}/lib/rosbag/play".format(distro)
    cmd = [play, bag, "--clock", "-r", rate, "-s", start]
    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(main())
