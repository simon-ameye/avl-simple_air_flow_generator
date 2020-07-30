# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:24:31 2020
@author: Simon Ameye AVL AST FRANCE
"""
#--------------------Libraries--------------------
import numpy as np
import tkinter as tk
from stl import mesh
from tkinter import filedialog
root = tk.Tk()
root.withdraw()


#--------------------User Parameters--------------------

# Data

file_path =  filedialog.askopenfilename(initialdir = "file_path",title = "Select file",filetypes = (("STL files","*.STL"),("all files","*.*")))
Gear_Mesh = mesh.Mesh.from_file(file_path)

Thickness_Points_Ratio = 0.05 #Ratio of furthest points from rotation axis to consider for the determination of gear thickness
Helix_angle = 0             /180*np.pi     #rad Default : 0
Wheel_Rotation_Speed = -226.1725                #RPS
Rotation_Direction = -1
Scaling_Of_Mesh = 1/1000
Domain_Scaling_To_Radius = 2

#Grid properties
Spacial_Discretization = 3   /1000          #m Recommended : Spacial_Discretization = 5

#Model properties
Dumping_oefficient = 1000                   #tuning Default : 1000


#Rescale mesh
Gear_Mesh[:,:] = Gear_Mesh[:,:]*Scaling_Of_Mesh

#Get important data from mesh
volume, Center_Of_Gravity, inertia = Gear_Mesh.get_mass_properties()

# Find rotation axis : If X and Y inearia are the same, Z is the revolution axis
# First, we isolate main inertias
Main_Axis_Inertias = np.array([[True,False,False],[False,True,False],[False,False,True]])
Inertias_XYZ = inertia[Main_Axis_Inertias]
#Then, we find the value that is the most different from the 2 others, ie the furthest from average
Revolution_Axis = np.argmax(abs(Inertias_XYZ-np.average(Inertias_XYZ))) # 0 for x, 1 for y, 2 for z

#If rot. axis = 0, we will check the radius accordinf to dir. 0 and 2. We isolate the other directions
Other_Dir_Than_Rot = np.array([0,1,2])[[0,1,2]!=Revolution_Axis]

#Find radius
#First, we isolate the first list of point of the mesh.
List_Of_Points = Gear_Mesh[:,0:3]
#Then, we look at all their distance to revolution axis
List_Of_All_Radius_From_Revol_Axis = np.sqrt(np.square(List_Of_Points[:,Other_Dir_Than_Rot[0]]-Center_Of_Gravity[Other_Dir_Than_Rot[0]])   +   np.square(List_Of_Points[:,Other_Dir_Than_Rot[1]]-Center_Of_Gravity[Other_Dir_Than_Rot[1]]))
#The radius is the max of all those distances
Gear_Radius = np.amax(List_Of_All_Radius_From_Revol_Axis)

#Find thickness limits of the gear
#We select the most furthest points from axis according to the Thickness_Points_Ratio
Points_Of_The_Edge_Of_Gear = List_Of_All_Radius_From_Revol_Axis>(Gear_Radius*(1-Thickness_Points_Ratio))
#Then, we look at their extreme positions
Gear_Extreme_Thickness_Pos = [np.amin(List_Of_Points[:,Revolution_Axis][Points_Of_The_Edge_Of_Gear]), np.amax(List_Of_Points[:,Revolution_Axis][Points_Of_The_Edge_Of_Gear])]

#Let's plot some data
Revolution_Axis_String = ['X','Y','Z'][Revolution_Axis]
print("Volume                                  = {0}".format(volume))
print("Position of the center of gravity (COG) = {0}".format(Center_Of_Gravity)+"m")
print("Inertia matrix at expressed at the COG  = {0}".format(inertia[0,:]))
print("                                          {0}".format(inertia[1,:]))
print("                                          {0}".format(inertia[2,:]))
print("Revolution axis = "+Revolution_Axis_String)
print("Gear radius = "+str(Gear_Radius)+"m")
print("Gear tooth limits according to "+Revolution_Axis_String+" = "+str(Gear_Extreme_Thickness_Pos)+"m")
print("Gear thickness = "+str(Gear_Extreme_Thickness_Pos[1]-Gear_Extreme_Thickness_Pos[0])+"m")


#Domain properties
Domain_Radius = Gear_Radius*Domain_Scaling_To_Radius
Gear_Thickness = Gear_Extreme_Thickness_Pos[1]-Gear_Extreme_Thickness_Pos[0]                                          #m Recommended : Domain_Radius = Wheel_Radius*2
Domain_thickness_limits = [Gear_Extreme_Thickness_Pos[0]-Gear_Thickness,Gear_Extreme_Thickness_Pos[1]+Gear_Thickness]            #m Recommended : Domain_thickness = (Domain_Radius-Wheel_Radius)*2+Wheel_thickness


#--------------------Code--------------------
#Disc sampling method
num_pts_in_disc = np.floor(np.pi*np.square(Domain_Radius/Spacial_Discretization))
indices = np.arange(0, num_pts_in_disc, dtype=float)
radius_list_init = Domain_Radius*np.sqrt(indices/num_pts_in_disc)
angle_list_init = np.pi * (1 + 5**0.5) * indices
axial_pos_list_init = np.arange(Domain_thickness_limits[0],Domain_thickness_limits[1],Spacial_Discretization)

#Lists initialization
angle_list = np.tile(angle_list_init,len(axial_pos_list_init))
radius_list = np.tile(radius_list_init,len(axial_pos_list_init))
axial_pos_list = np.repeat(axial_pos_list_init,len(angle_list_init))

#Result calculation
result = np.zeros([len(angle_list),6])
result[:,Revolution_Axis] = axial_pos_list
result[:,Other_Dir_Than_Rot[0]] = np.cos(angle_list)*radius_list+Center_Of_Gravity[Other_Dir_Than_Rot[0]]
result[:,Other_Dir_Than_Rot[1]] = np.sin(angle_list)*radius_list+Center_Of_Gravity[Other_Dir_Than_Rot[1]]

Points_Between_Gears_Axial_Extreme_Points = (axial_pos_list>=Gear_Extreme_Thickness_Pos[0])*(axial_pos_list<=Gear_Extreme_Thickness_Pos[1])
Axial_Distance_From_Gear_Extreme_Points = np.amin(np.abs([axial_pos_list-Gear_Extreme_Thickness_Pos[0],axial_pos_list-Gear_Extreme_Thickness_Pos[1]]),0)
axial_ratio = Points_Between_Gears_Axial_Extreme_Points*1            +         np.logical_not(Points_Between_Gears_Axial_Extreme_Points)*1/(1+np.square(Axial_Distance_From_Gear_Extreme_Points)*Dumping_oefficient)

radial_ratio = (radius_list<=Gear_Radius)*1 + (radius_list>Gear_Radius)*1/(1+np.square(radius_list-Gear_Radius)*Dumping_oefficient)
considered_radius = (radius_list<=Gear_Radius)*radius_list  +  (radius_list>Gear_Radius)*Gear_Radius
Air_Velocity = considered_radius*Wheel_Rotation_Speed*axial_ratio*radial_ratio
result[:,Other_Dir_Than_Rot[0]+3] = -np.sin(angle_list)*Air_Velocity
result[:,Other_Dir_Than_Rot[1]+3] = np.cos(angle_list)*Air_Velocity
result[:,Revolution_Axis+3] = Air_Velocity*np.sin(Helix_angle)*np.cos(Helix_angle)

#Write file
NAMES = np.array(['X','Y','Z','U','V','W']);
result_cat = np.vstack((NAMES,result))
file_path =  filedialog.asksaveasfilename(initialdir = "file_path",title = "Select file",filetypes = (("CSV files","*.CSV"),("all files","*.*")))
np.savetxt(file_path, result_cat, delimiter="," , fmt="%s")


