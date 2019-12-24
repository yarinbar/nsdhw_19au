#  Geographic Routing Simulator
## Background
Geographic Routing (GR), as its name, is an innovative method to utilize geographic information to complete routing task in Ad-Hoc network, especially Delay Tolerant Network (DTN) <sup>[1]</sup>. The main feature of DTN is that all nodes in this kind of network possess low computation ability, need to be very power-efficient (which means most of time they need to be in power-saved mode or even in sleep mode and weak up regularly) and the memory of them is too little to store traditional routing table, full network architecture, etc. 
Recently, with the great development of IoT, Ad-Hoc become a more and more popular topic in univesities or company’s development center. Many researches who study on network routing protocols or MAC layer architectures need a tool to complete the simulation task and validate their algorithm. However, based on my previous experience, there doesn’t exist a network simulation tool which is designed exactly for network that uses geographic information. At the same time, in modern IoT application, DTN plays a very important row in many actual scenario, such as outer-space networks, under-water networks, vehicular networks and sparse mobile networks which are featured in data exchanging happened in a long period and nodes are high mobile <sup>[2]</sup>. 
Thus, I consider to develop a special network simulator for GR which maybe pretty useful for many researcher in this field.
## Related work
In computer network research, network simulation is a technique whereby a software program models the behavior of a network by calculating the interaction between the different network entities (routers, switches, nodes, access points, links etc.).<sup>[3]</sup> There are several popular open-source or commercial network simulators. For example, OPNET, Tetcos and NS2/NS3 (except the last one, others are commercial software).
There’s also some simulators which are not develop for GR but can be modified to satisfy the basic function. For example, one-simulator, a opportunistic network routing simulator.
However, as I described in last sector, all of these popular simulator donot well support GR and need additional plug-in to execute simulation. The plug-in are not well compatible with software and cause bugs every time. I believe it’s very emerging to develop an open-source tool for related researchers. 
## Functions and TODO List

+ Network Nodes object
+ Nodes movement simulation
+ Log recording function
+ Implement several mainstream GR routing protocol
+ Iterator
+ ...

## Project Link (Github)
[GR simulator (Github)](http://example.com https://github.com/liuyuchen777/GR-simulator)
## Referece
[1]  C. Perkins, E. Belding-Royer, and S. Das, “Rfc 3561 ad hoc on-demand distance vector (aodv) routing,” United States, 2003.
[2] Shyam Kapadia, Bhaskar Krishnamachari and Lin Zhang, ‘Data Delivery in Delay Tolerant Networks: A Survey’, Cisco Systems Inc., San Jose, CA.
[3] Wikipedia contributors. (2019, September 8). Network simulation. In Wikipedia, The Free Encyclopedia. Retrieved 10:57, October 23, 2019.
[4] [The One simulator (an opportunistic network simulator)](http://example.com https://www.netlab.tkk.fi/tutkimus/dtn/theone/)
[5] [NS3 (The Network Simulator 3)](http://example.com https://www.nsnam.org)