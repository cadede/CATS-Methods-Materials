# CATS-Methods-Materials

Biologging tags are capable of capturing high-resolution data from cryptic species that spend a significant portion of their lives in remote or inaccessible environments, such as below the sea surface. The advancement of biologging technologies over the last ~30 years has lead to an explosion of research into the hidden lives of these elusive animals. Information about locomotion, feeding, migratory patterns, reproduction, and metabolism have all been determined using data obtained with biologging tags. CATS (customized animal tracking solutions) tags are one of the newest platforms for this type of research and integrate high-resolution video with a suite of inertial and environment sensors.

For the last 5 years, we as members of the Goldbogen lab at Stanford University's Hopkins Marine Station have been using CATS tags to study the biomechanics and behavior of large baleen whale species. In the course of this work, we have tested and streamlined a workflow and set of tools to integrate multiple streams of raw sensor data into a simple data packet that can be used for further analyses.

The tool kit works primarily with Matlab and includes tools to read/write, calibrate, process, visualize, and carry out statistical analysis of datasets. These tools are based on many of the tools designed at animaltags.org(link is external), but extend these ideas into a single tag processing workflow. The tools can work with a variety of tags (e.g. Acousonde, openTags, TDR-10) and are focused on CATS tags as models for video/data integration. The basic processing framework is also suitable for a variety of species, though our lab focuses primarily on large cetaceans.

This repository is designed to provide the entire CATS toolkit and a full tutorial (available in the Wiki section of this repository) for working with the CATS tag plaform. Our goal is for you to be able to turn raw data into a calibrated and properly oriented data packet that can be used for any of your subsequent analysis.

1) To start, download the CATS Tools folder and **add it to your Matlab file path with all folders and subfolders.**

2) If you have a brand new CATS tag, navigate to the Wiki and start with the first section ("Setting Up Tag And User Interface") to set up your tag. Please note that the CATS tag design changes regularly and may make some segments of the Wiki out-of-date. We will try to update the Wiki as often as we can.

3) If you have already set up your tag and successfully deployed it, navigate to the fourth section ("Importing CATS Data") and start working along from there.
