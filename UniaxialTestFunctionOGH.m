function [X,Y]=UniaxialTestFunctionOGH(deleteFile,sampleDim,sampleMesh, c, k1, k2, gamma, kappa, k, thickness_ini, BC, strain, ratio, rho_mat)

    %%%This is a test for The inverse uniaxial test using a multi layerd cube

    %%
%     clear; close all; clc;

    %% Plot settings
    fontSize=20;
    faceAlpha1=0.8;
    markerSize=40;
    markerSize2=35;
    lineWidth=3;
    cMap= 'jet'; %viridis(20); %colormap

    %% Control parameters

    % Path names

    savePath=fullfile('C:','Users','ahmad','Documents','Biomechanics','GibbonMessAround');

    % Defining file names
    febioFebFileNamePart='UniaxialModel_OGH_febio';
    febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
    febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
    febioXpltFileName=[febioFebFileNamePart,'.xplt']; %FEBio xplt file name
    febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
    febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
    febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
    febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting principal stress
    
    if deleteFile == 1
    delete([febioFebFileName]);
    delete([febioLogFileName]);
    delete([febioLogFileName_disp]);
    delete([febioLogFileName_force]);
    delete([febioLogFileName_stress]);
    delete([febioLogFileName_stress_prin]);
    delete([febioXpltFileName]);   
    end
    


    %Specifying dimensions and number of elements
    % cubeSize=10;
    % sampleWidth=cubeSize; %Width
    % sampleThickness=cubeSize; %Thickness
    % sampleHeight=cubeSize; %Height
    % pointSpacings=2*ones(1,3); %Desired point spacing between nodes
    % numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
    % numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
    % numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3

    sampleWidth= sampleDim(1); %1.5; %Width ie x axis
    sampleThickness=sampleDim(2); %1.5; %Thickness ie y axis
    sampleHeight=sampleDim(3); %1.5; %Height ie z axis

    numElementsWidth= sampleMesh(1); %10; %Number of elemens in dir 1
    numElementsThickness= sampleMesh(2); %10; %Number of elemens in dir 2
    numElementsHeight_TotalLayer= sampleMesh(3); %10; %Number of big layers in dir 3
    numElementsHeight_PerLayer= sampleMesh(4); %1; %Number of layers per big layer

    % Define applied displacement
    appliedStrain= strain; %Linear strain (Only used to compute applied stretch)
    loadingOption='tension'; % or 'tension'
    switch loadingOption
        case 'compression'
            stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
        case 'tension'
            stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
    end
    displacementMagnitude=(stretchLoad*sampleWidth)-sampleWidth; %The displacement magnitude

    % stretchLoad=0.7;
    % displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight;



    %Material parameter set

    %The inital guess parameters
%     E_ini = [1.174 3.522 2.348];
%     thickness_ini = [0 0.675 1.5]; %[0.15 0.6750 1.2750];

    %The interpolation to a smooth curve
    thicnkness_used = linspace((sampleHeight/numElementsHeight_TotalLayer)/2,sampleHeight-((sampleHeight/numElementsHeight_TotalLayer)/2),numElementsHeight_TotalLayer);
    
    c_layer = interp1(thickness_ini,c,thicnkness_used,'pchip');
    k1_layer = interp1(thickness_ini,k1,thicnkness_used,'pchip');
    k2_layer = interp1(thickness_ini,k2,thicnkness_used,'pchip');
    gamma_layer = interp1(thickness_ini,gamma,thicnkness_used,'pchip');
    kappa_layer = interp1(thickness_ini,kappa,thicnkness_used,'pchip');
    k_layer = interp1(thickness_ini,k,thicnkness_used,'pchip');
    
    
%     rho_mat = 1.12; %Material Density
        
    cFigure; hold on;
    plot(thickness_ini,c,'--k','lineWidth',lineWidth);
    plot(thicnkness_used,c_layer,'-o','lineWidth',lineWidth);
    xlabel('Thickness (m)');
    ylabel("c (Pa)");
    title("Material acros thickness",'FontSize',fontSize);
    legend("Initial","Used");

    view(2); axis tight;  grid on; grid minor; axis square; box on;
    set(gca,'FontSize',fontSize);
    drawnow;



    % FEA control settings
    numTimeSteps=10;%10; %Number of time steps desired
    max_refs=5; %25; %Max reforms
    max_ups=50; %0; %Set to zero to use full-Newton iterations
    opt_iter=6; %Optimum number of iterations
    max_retries=5; %Maximum number of retires
    dtmin=(1/numTimeSteps)/100; %Minimum time step size
    dtmax=1/numTimeSteps; %Maximum time step size
    runMode='internal';

    %% SIMULATE EXPERIMENTAL DATA

    % %Basic set
    % stress_cauchy_exp=1/1000*[-0.606636933451196;-0.594598753306976;-0.582704841004989;-0.571357405135258;-0.560202987257958;-0.549116632489736;-0.538518403222691;-0.528087294560408;-0.518193056737126;-0.508206114096577;-0.498701595140669;-0.489855637164223;-0.480813541456146;-0.472386398119889;-0.463619435755875;-0.455563887366101;-0.447492483369391;-0.439573886089611;-0.432050298442763;-0.424607647116797;-0.416804189884078;-0.410387298955262;-0.402977977822379;-0.396396657790034;-0.389210485373911;-0.383000553144204;-0.376675743693335;-0.370668858911072;-0.364731155035823;-0.358344772157269;-0.352790185960043;-0.346625957990168;-0.340956058045645;-0.335892515500584;-0.330212348100342;-0.325153422018813;-0.319890421672462;-0.315056500840712;-0.310859570288282;-0.305563240532117;-0.301114864342368;-0.295807178919732;-0.291944875824590;-0.287799721606394;-0.282704271932097;-0.279560319546267;-0.273953092186896;-0.271205596632553;-0.266019580975468;-0.261921529885230;-0.259473236771767;-0.254229845700605;-0.251227010966108;-0.246731599709182;-0.243347463269765;-0.240668206009318;-0.235904450179518;-0.233443491646300;-0.229240342796589;-0.226328455230997;-0.222574693739149;-0.219690552720043;-0.215908110296801;-0.213462994691799;-0.209402262394587;-0.206143135063048;-0.204259473767410;-0.200271046174199;-0.198497342254049;-0.194018107075590;-0.190682588685824;-0.190178278993820;-0.184939186637633;-0.184540226448861;-0.179325520197559;-0.177302998325867;-0.174896317893232;-0.170891038492450;-0.170506389072493;-0.165503062182587;-0.164964944739691;-0.160899776454826;-0.158388071874370;-0.156732253086585;-0.152865980799647;-0.151886036142296;-0.147064551962397;-0.146636586148680;-0.143247545748075;-0.139910407552933;-0.139643630040939;-0.135175245456319;-0.134411814767664;-0.131535143940800;-0.127943005303573;-0.127499404828055;-0.123718865018965;-0.123269655840332;-0.118450118919226;-0.117869603457104;-0.114259063948408;-0.111845005007273;-0.110782903827826;-0.106815200840467;-0.108112322079051;-0.103218831561054;-0.103859461792770;-0.100330051927225;-0.0988503888488038;-0.0984110683795259;-0.0920373613042230;-0.0944900398318279;-0.0908054642234128;-0.0873647791392896;-0.0857302637239363;-0.0832930518728098;-0.0811377337125286;-0.0801419455213994;-0.0773146108678843;-0.0750524119380378;-0.0737660915109812;-0.0711063097725948;-0.0689106003957611;-0.0662015603338655;-0.0637907034798034;-0.0622238776663924;-0.0587129121234732;-0.0590737570248270;-0.0542752113831988;-0.0539468997803651;-0.0504474583208646;-0.0479308792263506;-0.0474997497002284;-0.0422136232687380;-0.0419340474843669;-0.0383206523546593;-0.0353822402853126;-0.0342394575632298;-0.0296092241247699;-0.0290386117855990;-0.0252785740102147;-0.0211393477778685;-0.0210232271972257;-0.0149625128602809;-0.0150455267730763;-0.00925788965002460;-0.00559693887219605;-0.00235368730112040;0.00439939147625970;0.00280776088737496];
    % stretch_exp=[0.700330019000000;0.702340563275168;0.704351107550336;0.706361651825503;0.708372196100671;0.710382740375839;0.712393284651007;0.714403828926175;0.716414373201342;0.718424917476510;0.720435461751678;0.722446006026846;0.724456550302013;0.726467094577181;0.728477638852349;0.730488183127517;0.732498727402685;0.734509271677852;0.736519815953020;0.738530360228188;0.740540904503356;0.742551448778524;0.744561993053691;0.746572537328859;0.748583081604027;0.750593625879195;0.752604170154362;0.754614714429530;0.756625258704698;0.758635802979866;0.760646347255034;0.762656891530201;0.764667435805369;0.766677980080537;0.768688524355705;0.770699068630873;0.772709612906040;0.774720157181208;0.776730701456376;0.778741245731544;0.780751790006711;0.782762334281879;0.784772878557047;0.786783422832215;0.788793967107383;0.790804511382550;0.792815055657718;0.794825599932886;0.796836144208054;0.798846688483222;0.800857232758389;0.802867777033557;0.804878321308725;0.806888865583893;0.808899409859060;0.810909954134228;0.812920498409396;0.814931042684564;0.816941586959732;0.818952131234899;0.820962675510067;0.822973219785235;0.824983764060403;0.826994308335570;0.829004852610738;0.831015396885906;0.833025941161074;0.835036485436242;0.837047029711409;0.839057573986577;0.841068118261745;0.843078662536913;0.845089206812081;0.847099751087248;0.849110295362416;0.851120839637584;0.853131383912752;0.855141928187920;0.857152472463087;0.859163016738255;0.861173561013423;0.863184105288591;0.865194649563758;0.867205193838926;0.869215738114094;0.871226282389262;0.873236826664430;0.875247370939597;0.877257915214765;0.879268459489933;0.881279003765101;0.883289548040269;0.885300092315436;0.887310636590604;0.889321180865772;0.891331725140940;0.893342269416107;0.895352813691275;0.897363357966443;0.899373902241611;0.901384446516779;0.903394990791946;0.905405535067114;0.907416079342282;0.909426623617450;0.911437167892617;0.913447712167785;0.915458256442953;0.917468800718121;0.919479344993289;0.921489889268456;0.923500433543624;0.925510977818792;0.927521522093960;0.929532066369128;0.931542610644295;0.933553154919463;0.935563699194631;0.937574243469799;0.939584787744967;0.941595332020134;0.943605876295302;0.945616420570470;0.947626964845638;0.949637509120805;0.951648053395973;0.953658597671141;0.955669141946309;0.957679686221477;0.959690230496644;0.961700774771812;0.963711319046980;0.965721863322148;0.967732407597316;0.969742951872483;0.971753496147651;0.973764040422819;0.975774584697987;0.977785128973154;0.979795673248322;0.981806217523490;0.983816761798658;0.985827306073826;0.987837850348993;0.989848394624161;0.991858938899329;0.993869483174497;0.995880027449664;0.997890571724832;0.999901116000000];
    %
    % %Interpolate to higher sampling
    % n=100;
    % stretch_exp_n=linspace(1,stretchLoad,n);
    % stress_cauchy_exp_n = interp1(stretch_exp,stress_cauchy_exp,stretch_exp_n,'pchip');
    %
    % %Override variables
    % stress_cauchy_exp=stress_cauchy_exp_n;
    % stretch_exp=stretch_exp_n;
    %
    % %Add noise
    % stdNoise=0.01; %Standard deviation in units of stress
    % stress_cauchy_exp_n=stress_cauchy_exp_n+stdNoise.*randn(size(stress_cauchy_exp_n));



    %% Creating model geometry and mesh
    % A box is created with tri-linear hexahedral (hex8) elements using the
    % |hexMeshBox| function. The function offers the boundary faces with
    % seperate labels for the top, bottom, left, right, front, and back sides.
    % As such these can be used to define boundary conditions on the exterior.

    % Create a box with hexahedral elements
    cubeDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
    cubeElementNumbers=[numElementsWidth numElementsThickness numElementsHeight_TotalLayer numElementsHeight_PerLayer]; %Number of elements
    outputStructType=2; %A structure compatible with mesh view
    [meshStruct]=hexMeshBoxPerLayer(cubeDimensions,cubeElementNumbers,outputStructType);

    %Access elements, nodes, and faces from the structure
    E=meshStruct.elements; %The elements
    E_Layer = meshStruct.elements_layer;
    V=meshStruct.nodes; %The nodes (vertices)
    Fb=meshStruct.facesBoundary; %The boundary faces
    Cb=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
    elementMaterialIndices=ones(size(E,1),1); %Element material indices


    c_temp = [];
    for k_loop = 1:numElementsHeight_TotalLayer

        c_temp = [c_temp; ones(size(E_Layer,1),1)*c_layer(k_loop)];
    end
    meshStruct.elementData= c_temp;
    %%
    % Plotting model boundary surfaces and a cut view

    hFig=cFigure;

    subplot(1,2,1); hold on;
    title('Model boundary surfaces and labels','FontSize',fontSize);
    gpatch(Fb,V,Cb,'k',faceAlpha1);
    colormap(gjet(6)); icolorbar;
    axisGeom(gca,fontSize);

    hs=subplot(1,2,2); hold on;
    title('Cut view of solid mesh','FontSize',fontSize);
    optionStruct.hFig=[hFig hs];
    optionStruct.numSLiceSteps=numElementsHeight_TotalLayer/numElementsHeight_PerLayer;
    optionStruct.cutDir=3;
    optionStruct.cutSide=-1;
    meshView(meshStruct,optionStruct);
    colormap(gca);
    axisGeom(gca,fontSize);

    drawnow;

    %% Defining the boundary conditions
    % The visualization of the model boundary shows colors for each side of the
    % cube. These labels can be used to define boundary conditions.

    %Define supported node sets
    logicFace=Cb==1; %Logic for current face set
    Fr=Fb(logicFace,:); %The current face set
    bcSupportList_X=unique(Fr(:)); %Node set part of selected face

    logicFace=Cb==3; %Logic for current face set
    Fr=Fb(logicFace,:); %The current face set
    bcSupportList_Y=unique(Fr(:)); %Node set part of selected face

    logicFace=Cb==5; %Logic for current face set
    Fr=Fb(logicFace,:); %The current face set
    bcSupportList_Z=unique(Fr(:)); %Node set part of selected face
      

    %Prescribed displacement nodes
    logicPrescribe=Cb==2; %Logic for current face set
    Fr=Fb(logicPrescribe,:); %The current face set
    bcPrescribeList_X=unique(Fr(:)); %Node set part of selected face
    
    logicPrescribe=Cb==4; %Logic for current face set
    Fr=Fb(logicPrescribe,:); %The current face set
    bcPrescribeList_Y=unique(Fr(:)); %Node set part of selected face

    %%
    % Visualizing boundary conditions. Markers plotted on the semi-transparent
    % model denote the nodes in the various boundary condition lists.

    hf=cFigure;
    title('Boundary conditions','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hold on;

    gpatch(Fb,V,'kw','k',0.5);

    hl(1)=plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize);
    hl(2)=plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize);
    hl(3)=plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize);
    hl(4)=plotV(V(bcPrescribeList_X,:),'k.','MarkerSize',markerSize);
    hl(5)=plotV(V(bcPrescribeList_Y,:),'m.','MarkerSize',markerSize);

    legend(hl,{'BC_X support','BC_Y support','BC_Z support','Prescribe_X','Prescribe_Y'});
    % legend(hl,{'BC support','BC Prescribed'});

    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;

    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.

    %Get a template with default settings
    [febio_spec]=febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version='3.0';

    %Module section
    febio_spec.Module.ATTR.type='solid';

    %Control section
    febio_spec.Control.analysis='DYNAMIC'; %'STATIC';
    febio_spec.Control.time_steps=numTimeSteps;
    febio_spec.Control.step_size=1/numTimeSteps;
    febio_spec.Control.solver.max_refs=max_refs;
    febio_spec.Control.solver.max_ups=max_ups;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;
    

   febio_spec.Control.solver.etol=0.01;    
   febio_spec.Control.solver.dtol=0.001;    
   febio_spec.Control.solver.rtol=0.001;
   febio_spec.Control.solver.lstol=0.9;
   febio_spec.Control.solver.min_residual=1e-20;
   febio_spec.Control.solver.qnmethod='BROYDEN'; 
   febio_spec.Control.solver.rhoi=0;
   febio_spec.Control.solver.symmetric_stiffness=0;
   
    febio_spec.Globals.Constants.T = 0;
    febio_spec.Globals.Constants.R = 0;
    febio_spec.Globals.Constants.Fc = 0;
    
    %Material section
    
% <material id="2" type="Holzapfel-Gasser-Ogden">
%   <c>7.64</c>
%   <k1>996.6</k1>
%   <k2>524.6</k2>
%   <gamma>49.98</gamma>
%   <kappa>0.226</kappa>
%   <k>1e5</k>
% </material>

    for k_loop = 1:numElementsHeight_TotalLayer
        febio_spec.Material.material{k_loop}.ATTR.name= ['Material',num2str(k_loop)]; %char(materialName(k));
        febio_spec.Material.material{k_loop}.ATTR.type='Holzapfel-Gasser-Ogden';
        febio_spec.Material.material{k_loop}.ATTR.id=k_loop;
        febio_spec.Material.material{k_loop}.c =c_layer(k_loop);
        febio_spec.Material.material{k_loop}.k1 =k1_layer(k_loop);
        febio_spec.Material.material{k_loop}.k2 =k2_layer(k_loop);
        febio_spec.Material.material{k_loop}.gamma =gamma_layer(k_loop);
        febio_spec.Material.material{k_loop}.kappa =kappa_layer(k_loop);
        febio_spec.Material.material{k_loop}.k =k_layer(k_loop);

        febio_spec.Material.material{k_loop}.density =rho_mat;
    end

    % Mesh section
    % -> Nodes
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

    E_id_start = 1;
    E_id_end = 0;
    T = [];

    % -> Elements
    % partName={'Part1','Part2','Part3','Part4','Part5'};
    for k_loop = 1:numElementsHeight_TotalLayer
        febio_spec.Mesh.Elements{k_loop}.ATTR.name= ['Part',num2str(k_loop)];%char(partName(k)); %Name of this part
        febio_spec.Mesh.Elements{k_loop}.ATTR.type='hex8'; %Element type

        %     febio_spec.Mesh.Elements{k}.elem.ATTR.id =(size(E_Layer,1)*(k-1))+(1:1:size(E_Layer,1))'; %Element id's
        E_id_end = E_id_end + size(E_Layer(:,:,k_loop),1);
        febio_spec.Mesh.Elements{k_loop}.elem.ATTR.id =(E_id_start:1:E_id_end)'; %Element id's
        E_id_start = E_id_end+1;

        febio_spec.Mesh.Elements{k_loop}.elem.VAL = E_Layer(:,:,k_loop); %The element matrix

        %     T(:,k) =(E_id_start:1:E_id_end)';
    end

    % -> NodeSets
    nodeSetName1='bcSupportList_X';
    nodeSetName2='bcSupportList_Y';
    nodeSetName3='bcSupportList_Z';
    nodeSetName4='bcPrescribeList_X';
    nodeSetName5='bcPrescribeList_Y';

    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList_X(:);

    febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
    febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcSupportList_Y(:);
    %
    febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
    febio_spec.Mesh.NodeSet{3}.node.ATTR.id=bcSupportList_Z(:);

    febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
    febio_spec.Mesh.NodeSet{4}.node.ATTR.id=bcPrescribeList_X(:);
    
    febio_spec.Mesh.NodeSet{5}.ATTR.name=nodeSetName5;
    febio_spec.Mesh.NodeSet{5}.node.ATTR.id=bcPrescribeList_Y(:);

    %MeshDomains section
    for k_loop = 1:numElementsHeight_TotalLayer
        febio_spec.MeshDomains.SolidDomain{k_loop}.ATTR.name = ['Part',num2str(k_loop)];%char(partName(k));
        febio_spec.MeshDomains.SolidDomain{k_loop}.ATTR.mat = ['Material',num2str(k_loop)];%char(materialName(k));
    end

    %Boundary condition section
    if BC ==1
        % -> Fix boundary conditions
        febio_spec.Boundary.bc{1}.ATTR.type='fix';
        febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
        febio_spec.Boundary.bc{1}.dofs='x';

        febio_spec.Boundary.bc{2}.ATTR.type='fix';
        febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
        febio_spec.Boundary.bc{2}.dofs='y';

        febio_spec.Boundary.bc{3}.ATTR.type='fix';
        febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
        febio_spec.Boundary.bc{3}.dofs='z';

        febio_spec.Boundary.bc{4}.ATTR.type='prescribe';
        febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
        febio_spec.Boundary.bc{4}.dof='x';
        febio_spec.Boundary.bc{4}.scale.ATTR.lc=1;
        febio_spec.Boundary.bc{4}.scale.VAL=displacementMagnitude;
        febio_spec.Boundary.bc{4}.relative=0;
        
        if ratio ~= 0
            febio_spec.Boundary.bc{5}.ATTR.type='prescribe';
            febio_spec.Boundary.bc{5}.ATTR.node_set=nodeSetName5;
            febio_spec.Boundary.bc{5}.dof='y';
            febio_spec.Boundary.bc{5}.scale.ATTR.lc=1;
            febio_spec.Boundary.bc{5}.scale.VAL=displacementMagnitude*ratio;
            febio_spec.Boundary.bc{5}.relative=0;
        end
    
    else
        
            % -> Fix boundary conditions
    febio_spec.Boundary.bc{1}.ATTR.type='fix';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.dofs='x';

    febio_spec.Boundary.bc{2}.ATTR.type='fix';
    febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{2}.dofs='y';

    febio_spec.Boundary.bc{3}.ATTR.type='fix';
    febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{3}.dofs='z';

    febio_spec.Boundary.bc{4}.ATTR.type='prescribe';
    febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
    febio_spec.Boundary.bc{4}.dof='x';
    febio_spec.Boundary.bc{4}.scale.ATTR.lc=1;
    febio_spec.Boundary.bc{4}.scale.VAL=displacementMagnitude;
    febio_spec.Boundary.bc{4}.relative=0;
      
    
    end


    %LoadData section
    % -> load_controller
    febio_spec.LoadData.load_controller{1}.ATTR.id=1;
    febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

    %Output section
    % -> log file
    febio_spec.Output.logfile.ATTR.file=febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='sx;sy;sz';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
    
    febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stress_prin;
    febio_spec.Output.logfile.element_data{2}.ATTR.data='s1;s2;s3';
    febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';


    %%%%%%%%%
    %%
    %%

    %% Quick viewing of the FEBio input file structure
    % The |febView| function can be used to view the xml structure in a MATLAB
    % figure window.

    %%
    % |febView(febio_spec); %Viewing the febio file|

    %% Exporting the FEBio input file
    % Exporting the febio_spec structure to an FEBio input file is done using
    % the |febioStruct2xml| function.

    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    %system(['gedit ',febioFebFileName,' &']);

    %% Running the FEBio analysis
    % To run the analysis defined by the created FEBio input file the
    % |runMonitorFEBio| function is used. The input for this function is a
    % structure defining job settings e.g. the FEBio input file name. The
    % optional output runFlag informs the user if the analysis was run
    % succesfully.

    febioAnalysis.run_filename=febioFebFileName; %The input file name
    febioAnalysis.run_logname=febioLogFileName; %The name for the log file
    febioAnalysis.disp_on=1; %Display information on the command window
    febioAnalysis.runMode=runMode;
    febioAnalysis.maxLogCheckTime=10; %Max log file checking time

    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    %% Import FEBio results

    if runFlag==1 %i.e. a succesful run

        %%

        % Importing nodal displacements from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);

        %Access data
        N_disp_mat=dataStruct.data; %Displacement
        timeVec=dataStruct.time; %Time

        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure
        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('Displacement magnitude [mm]','Interpreter','Latex')
        hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1,2); %Add graphics object to animate
        hp.Marker='.';
        hp.MarkerSize=markerSize2;
        hp.FaceColor='interp';
        gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0 max(DN_magnitude)]); caxis manual;
        axis(axisLim(V_DEF)); %Set axis limits statically
        view(140,30);
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments
            DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        %%
        % Importing element stress from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);

        %Access data
        E_stress_mat_temp=dataStruct.data;
        
        E_stress_mat_x = E_stress_mat_temp(:,1,:);
        E_stress_mat_y = E_stress_mat_temp(:,2,:);
        E_stress_mat_z = E_stress_mat_temp(:,3,:);

        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

%         [CV]=faceToVertexMeasure(E,V,E_stress_mat_x(:,:,end));

        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),1,1);
        E_stress_prin_mat=dataStruct.data;
        [CV]=faceToVertexMeasure(E,V,sqrt(sum(E_stress_prin_mat(:,:,end).^2,2)));
        

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy

        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('$\sigma_{Von Mises}$ [MPa]','Interpreter','Latex')
        hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1,2); %Add graphics object to animate
        hp.Marker='.';
        hp.MarkerSize=markerSize2;
        hp.FaceColor='interp';
        gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([min(E_stress_mat_x(:)) max(E_stress_mat_x(:))]);
        axis(axisLim(V_DEF)); %Set axis limits statically
        view(140,30);
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments

            [CV]=faceToVertexMeasure(E,V,E_stress_mat_x(:,:,qt));

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        %% Useless for now
        % Calculate metrics to visualize stretch-stress curve
        %
        %     DZ_set=N_disp_mat(bcPrescribeList,end,:); %Z displacements of the prescribed set
        %     DZ_set=mean(DZ_set,1); %Calculate mean Z displacements across nodes
        %     stretch_sim=(DZ_set(:)+sampleHeight)./sampleHeight; %Derive stretch
        %     stress_cauchy_sim=mean(squeeze(E_stress_mat(:,end,:)),1)';
        %
        %     % Visualize stress-stretch curve
        %
        %     cFigure; hold on;
        %     title('Uniaxial stress-stretch Z curve','FontSize',fontSize);
        %     xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
        %     ylabel('$\sigma_{zz}$ [MPa]','FontSize',fontSize,'Interpreter','Latex');
        %
        %     plot(stretch_sim(:),stress_cauchy_sim(:),'r-','lineWidth',lineWidth);
        %
        %     view(2); axis tight;  grid on; axis square; box on;
        %     set(gca,'FontSize',fontSize);
        %     drawnow;

        %% For the X direction New I've made now
        % Calculate metrics to visualize stretch-stress curve

        DX_set=N_disp_mat(bcPrescribeList_X,1,:); %X displacements of the prescribed set
        DX_set=mean(DX_set,1); %Calculate mean X displacements across nodes
        stretch_sim=(DX_set(:)+sampleWidth)./sampleWidth; %Derive stretch
        
        
        stress_cauchy_sim=mean(squeeze(E_stress_mat_x(:,end,:)),1)';
        
        %1000 1 11
        
%         for m = 1:size(E_stress_mat,3)
%            
%             stress_cauchy_sim(m) = mean(E_stress_mat(:,1,m));
%             
%         end
        
        
        cFigure; hold on;
        title('Uniaxial stress-strech X curve','FontSize',fontSize);
        xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
        ylabel('$\sigma$ [Pa]','FontSize',fontSize,'Interpreter','Latex');
%         ylabel('$\sigma_{xx}$ [MPa]','FontSize',fontSize,'Interpreter','Latex');

        plot(stretch_sim(:),stress_cauchy_sim(:),'r-','lineWidth',lineWidth);

        view(2); axis tight;  grid on; axis square; box on;
        set(gca,'FontSize',fontSize);
        drawnow;

%%
        % For the Y direction New I've made now
        % Calculate metrics to visualize stretch-stress curve

        DX_set=N_disp_mat(bcPrescribeList_X,1,:); %X displacements of the prescribed set
        DX_set=mean(DX_set,1); %Calculate mean X displacements across nodes
        stretch_sim=(DX_set(:)+sampleWidth)./sampleWidth; %Derive stretch
        
        
        stress_cauchy_sim=mean(squeeze(E_stress_mat_y(:,end,:)),1)';
        
        
        
%         cFigure; 
        hold on;
%         title('Uniaxial stress-strech Y curve','FontSize',fontSize);
        xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
%         ylabel('$\sigma_{yy}$ [MPa]','FontSize',fontSize,'Interpreter','Latex');

        plot(stretch_sim(:),stress_cauchy_sim(:),'g-','lineWidth',lineWidth);

        view(2); axis tight;  grid on; axis square; box on;
        set(gca,'FontSize',fontSize);
        drawnow;

        
        
        
%%
        % For the Z direction New I've made now
        % Calculate metrics to visualize stretch-stress curve

        DX_set=N_disp_mat(bcPrescribeList_X,1,:); %X displacements of the prescribed set
        DX_set=mean(DX_set,1); %Calculate mean X displacements across nodes
        stretch_sim=(DX_set(:)+sampleWidth)./sampleWidth; %Derive stretch
        
        
        stress_cauchy_sim=mean(squeeze(E_stress_mat_z(:,end,:)),1)';
        
        
        
%         cFigure; 
        hold on;
%         title('Uniaxial stress-strech Z curve','FontSize',fontSize);
        xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
%         ylabel('$\sigma_{zz}$ [MPa]','FontSize',fontSize,'Interpreter','Latex');

        plot(stretch_sim(:),stress_cauchy_sim(:),'b-','lineWidth',lineWidth);

        view(2); axis tight;  grid on; axis square; box on;
        set(gca,'FontSize',fontSize);
        
        
        drawnow;

        
%%
        % Importing element principal stresses from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),1,1);

        %Access data
        E_stress_prin_mat=dataStruct.data;
        
        P1t = E_stress_prin_mat(:,1,:);
        P2t = E_stress_prin_mat(:,2,:);
        P3t = E_stress_prin_mat(:,3,:);
        
        DX_set=N_disp_mat(bcPrescribeList_X,1,:); %X displacements of the prescribed set
        DX_set=mean(DX_set,1); %Calculate mean X displacements across nodes
        stretch_sim=(DX_set(:)+sampleWidth)./sampleWidth; %Derive stretch
        
        P1=mean(squeeze(P1t(:,end,:)),1)';
        P2=mean(squeeze(P2t(:,end,:)),1)';
        P3=mean(squeeze(P3t(:,end,:)),1)';
        
        PVM = sqrt( 0.5 *((P2-P3).^2 + (P3-P1).^2 + (P1-P2).^2));
        
%         cFigure; 
        hold on;
%         title('Vom Mises','FontSize',fontSize);
        xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
%         ylabel('$\sigma_{VM}$ [MPa]','FontSize',fontSize,'Interpreter','Latex');

        plot(stretch_sim(:),PVM(:),'k-','lineWidth',lineWidth);

        view(2); axis tight;  grid on; axis square; box on;
        set(gca,'FontSize',fontSize);
        
        legend({'$\sigma_X$','$\sigma_Y$','$\sigma_Z$','$\sigma_{VonMises}$'},'location','northwest','Interpreter','Latex');
        
        drawnow;
        
        
%%
        

%         %Compute pressure
%         P = squeeze(-1/3*mean(sum(E_stress_prin_mat,2),1));
% 
%         %%
%         % Visualize pressure-stretch curve
% 
%         cFigure; hold on;
%         title('Pressure-stretch curve','FontSize',fontSize);
%         xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex');
%         ylabel('$p$ [MPa]','FontSize',fontSize,'Interpreter','Latex');
% 
%         plot(stretch_sim(:),P(:),'r-','lineWidth',lineWidth);
% 
%         view(2); axis tight;  grid on; axis square; box on;
%         set(gca,'FontSize',fontSize);
%         drawnow;

    end



    %%
    %
    % <<gibbVerySmall.gif>>
    %
    % _*GIBBON*_
    % <www.gibboncode.org>
    %
    % _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>

    %%
    % _*GIBBON footer text*_
    %
    % License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
    %
    % GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
    % image segmentation, image-based modeling, meshing, and finite element
    % analysis.
    %
    % Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
    %
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
        Y = mean(squeeze(E_stress_mat_x(:,end,:)),1)';
        X =(DX_set(:)+sampleWidth)./sampleWidth;

end