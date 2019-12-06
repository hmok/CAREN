function [BodyPosGlob, BodyPosLoc, BodyCOMLoc, BodyCOMGlob, Masses, BWLoad, BodyNames, COMLoc] = Get_BodyProps_HMedit(modelpath)
%% Example: this is a modified verson of Get_BodyProps.m file
% [BodyPosGlob, BodyPosLoc, BodyCOMLoc, BodyCOMGlob, Masses, BWLoad, BodyNames, COMLoc] = Get_BodyProps_HMedit('10004_BL_SizeScaled.osim')
% example from prevoius version:
% Origanly written by Alex at Bouxsein Lab then I (Hossein M) made changes
% [BodyPosGlob, BodyPosLoc, BodyCOMLoc, BodyCOMGlob, Masses, BWLoad, BodyNames] = Get_BodyProps('10004_BL_SizeScaled.osim')

%% 
    import org.opensim.modeling.*

    myModel = Model(modelpath);
    state = myModel.initSystem;
    myModel.equilibrateMuscles(state);
    
    myBodies = myModel.getBodySet();
    numBodies = myBodies.getSize;
    BodyNames = cell(numBodies,1);
    
    for i=1:numBodies

        BodyNames{i} = char(myBodies.get(i-1).getName); %Array of Body names
        Masses(i,1) = myBodies.get(i-1).getMass();

        localPos = ArrayDouble.createVec3([0.0,0.0,0.0]);
        globalPos = ArrayDouble.createVec3([0.0,0.0,0.0]);
        myModel.getSimbodyEngine().getPosition(state, myBodies.get(i-1),localPos,globalPos);
        BodyPosGlob.(BodyNames{i}) = str2num(char(ArrayDouble.getValuesFromVec3(globalPos))); %Body location in global frame
        BodyPosLoc.(BodyNames{i}) = str2num(char(ArrayDouble.getValuesFromVec3(localPos)));   %Body location in local frame
        clear localPos globalPos

        localPos = ArrayDouble.createVec3([0,0,0]);
        globalPos = ArrayDouble.createVec3([0.0,0.0,0.0]);
%         myBodies.get(i-1).getMassCenter(localPos); 
        myBodies.get(i-1).getMassCenter();
        BodyCOMLoc.(BodyNames{i})= str2num(char(ArrayDouble.getValuesFromVec3(localPos)));   %Body COM location in local frame
        myModel.getSimbodyEngine().transformPosition(state, myBodies.get(i-1),localPos,myModel.getGround, globalPos);
        BodyCOMGlob.(BodyNames{i})= str2num(char(ArrayDouble.getValuesFromVec3(globalPos))); %Body COM location in global frame
        clear localPos globalPos

    end
    
    %load on vert due to body weight, T1 to L5
%     These are wsrong make changes later I changed to run them only 
    for i = 1:17 %L5 to T1
        BWLoad(i,1) = sum(Masses(i+5:end,1)); %BWLoad(i,1) = sum(Masses(i+5:22,1)); 
    end
    
    ArmMass = sum(Masses(:,1))*2;%ArmMass = sum(Masses(50:53,1))*2; 
    BWLoad = BWLoad + ArmMass;

    BWLoad(:,1) = BWLoad(end:-1:1,1); %re-orderBWLoad(:,1) = BWLoad(17:-1:1,1); %re-order
    BWLoad = BWLoad'*9.81;

    %% This is what Hossein added to get the location of center of mass for all upper body
    a= myModel.calcMassCenterPosition(state);
    COMLoc = str2num(char(ArrayDouble.getValuesFromVec3(a)));
    clear localPosCOM state
    
    
end