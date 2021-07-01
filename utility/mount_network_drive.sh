#!/bin/bash
mkdir -p [path_on_your_server_where_you_want_to_mount]
mount -t cifs -o vers=3.0,username=[username] [path_to_your_network_share_drive] [path_on_your_server_where_you_want_to_mount]