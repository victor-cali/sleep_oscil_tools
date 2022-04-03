classdef ThresholdReviewer_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        Panel              matlab.ui.container.Panel
        GridLayout2        matlab.ui.container.GridLayout
        GridLayout4        matlab.ui.container.GridLayout
        SaveButton         matlab.ui.control.Button
        RippleCheckBox     matlab.ui.control.CheckBox
        SharpwaveCheckBox  matlab.ui.control.CheckBox
        GridLayout8        matlab.ui.container.GridLayout
        RThresholdLabel    matlab.ui.control.Label
        RSpinner           matlab.ui.control.Spinner
        GridLayout7        matlab.ui.container.GridLayout
        SWThresholdLabel   matlab.ui.control.Label
        SWSpinner          matlab.ui.control.Spinner
        GridLayout         matlab.ui.container.GridLayout
        GridLayout3        matlab.ui.container.GridLayout
        EventsSlider       matlab.ui.control.Slider
        NextButton         matlab.ui.control.Button
        LastButton         matlab.ui.control.Button
        UIAxes             matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        x;
        bel;
        pyr;
        len;
        animal;
        filt_bel;
        filt_pyr;
        fn = 600;
        detections;
        window = 3601;
        detections_num;
        true_positive_ripples = [];
        true_positive_sharpwaves = [];
    end
    
    methods (Access = private)
        function update(app)
            half_wind = fix(app.window/2);
            % Get event number
            det_num = app.EventsSlider.Value;
            % start
            det_s = (app.detections.Peak(det_num) - half_wind);
            % end
            det_e = (app.detections.Peak(det_num) + half_wind);
            
            % Update saved or unsaved state
            if ismember(det_num, app.true_positive_ripples)
                app.RippleCheckBox.Value = 1;
            else
                app.RippleCheckBox.Value = 0;
            end
            if ismember(det_num, app.true_positive_sharpwaves)
                app.SharpwaveCheckBox.Value = 1;
            else
                app.SharpwaveCheckBox.Value = 0;
            end
            
            bel_sig = app.bel;
            space = max(bel_sig(det_s:det_e)) - min(bel_sig(det_s:det_e));
            pyr_sig = app.pyr + space;
            space = max(pyr_sig(det_s:det_e)) - min(bel_sig(det_s:det_e));
            filt_bel_sig = app.filt_bel + space;
            space = max(filt_bel_sig(det_s:det_e)) - min(bel_sig(det_s:det_e));
            filt_pyr_sig = 5.*app.filt_pyr + space; 
            
            y_max = max(filt_pyr_sig(det_s:det_e)) + 250;
            y_min = min(bel_sig(det_s:det_e)) - 250;

            title(app.UIAxes, strcat("Event ", int2str(det_num)))
            cla(app.UIAxes)
            plot(app.UIAxes, app.x, filt_pyr_sig, 'Color', 'blue')
            hold(app.UIAxes,'on')
            plot(app.UIAxes, app.x, filt_bel_sig, 'Color', 'black')
            hold(app.UIAxes,'on')
            plot(app.UIAxes, app.x, pyr_sig, 'Color', 'blue')
            hold(app.UIAxes,'on')
            plot(app.UIAxes, app.x, bel_sig, 'Color', 'black')
            
            app.UIAxes.YLim = [y_min y_max];
            app.UIAxes.XLim = [det_s/app.fn, det_e/app.fn];
            
            hold(app.UIAxes,'on')
            txt = app.detections.Type(det_num);
            text(app.UIAxes, det_e-1, y_max+250, txt)

        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            [file,path] = uigetfile('',...
                          'Select below and pyramidal signals', ... 
                          'MultiSelect', 'on');
            
            if length(file) == 2
                
                app.bel = importdata(fullfile(path, file{1}));
                app.pyr = importdata(fullfile(path, file{2}));
                app.len = min([length(app.bel) length(app.pyr)]);
                app.bel = app.bel(1:app.len);
                app.pyr = app.pyr(1:app.len);
                app.x = (1:app.len)/app.fn;
                
                [c,d] = butter(3, [1/300 20/300]);
                [e,f] = butter(3, [90/300 200/300]);
                app.filt_bel = filtfilt(c,d, app.bel);
                app.filt_pyr = filtfilt(e,f, app.pyr);
                
                [file,path] = uigetfile('','Select detections');
                
                if file ~= 0
                    
                    app.animal = strcat('_rat_', file(8:end-4));
                
                    filename = fullfile(path,file);
                    
                    variable = 'grouped_oscil_table';
                    S = load(filename, variable);
                    app.detections = S.(variable);
                    
                    app.detections_num = height(app.detections);
                    app.EventsSlider.Limits = [1 app.detections_num];

                    app.SaveButton.Enable = 'on';
                    app.SharpwaveCheckBox.Enable = 'on';
                    app.RippleCheckBox.Enable = 'on';
                    app.EventsSlider.Enable = 'on';
                    app.NextButton.Enable = 'on';
                    app.LastButton.Enable = 'on';
                    
                    app.update()
                    
                else
                    delete(app);
                end
            else
                delete(app);
            end
        end

        % Button pushed function: NextButton
        function NextButtonPushed(app, event)
            if app.EventsSlider.Value < app.detections_num
                app.EventsSlider.Value = app.EventsSlider.Value + 1;
                app.update()
            end
        end

        % Button pushed function: LastButton
        function LastButtonPushed(app, event)
            if app.EventsSlider.Value > 1
                app.EventsSlider.Value = app.EventsSlider.Value - 1;
                app.update()
            end
        end

        % Value changed function: EventsSlider
        function EventsSliderValueChanged(app, event)
            value = round(app.EventsSlider.Value);
            app.EventsSlider.Value = value;
            app.update()
        end

        % Value changed function: RippleCheckBox
        function RippleCheckBoxValueChanged(app, event)
            value = app.RippleCheckBox.Value;
            if value
                app.true_positive_ripples = [app.true_positive_ripples app.EventsSlider.Value];
            else
                app.true_positive_ripples(app.true_positive_ripples == app.EventsSlider.Value) = [];
            end
            app.update()
        end

        % Value changed function: SharpwaveCheckBox
        function SharpwaveCheckBoxValueChanged(app, event)
            value = app.SharpwaveCheckBox.Value;
            if value
                app.true_positive_sharpwaves = [app.true_positive_sharpwaves app.EventsSlider.Value];
            else
                app.true_positive_sharpwaves(app.true_positive_sharpwaves == app.EventsSlider.Value) = [];
            end
            app.update()
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            %ripples = app.true_positive_ripples;
            %sharpwaves = app.true_positive_sharpwaves;
            vars = {'app.true_positive_ripples', ...
                    'app.true_positive_sharpwaves'};
            file = strcat('threshold_review_', ...
                   erase(date,'-'), app.rat);
            uisave(vars, file);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'3x', '1x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Event')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 1;

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Layout.Row = 2;
            app.Panel.Layout.Column = 1;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel);
            app.GridLayout2.ColumnWidth = {'1x'};
            app.GridLayout2.RowHeight = {'1.5x', '1x'};

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout2);
            app.GridLayout3.ColumnWidth = {'1x', '8x', '1x'};
            app.GridLayout3.RowHeight = {'1x'};
            app.GridLayout3.Padding = [0 0 0 0];
            app.GridLayout3.Layout.Row = 2;
            app.GridLayout3.Layout.Column = 1;

            % Create LastButton
            app.LastButton = uibutton(app.GridLayout3, 'push');
            app.LastButton.ButtonPushedFcn = createCallbackFcn(app, @LastButtonPushed, true);
            app.LastButton.Enable = 'off';
            app.LastButton.Layout.Row = 1;
            app.LastButton.Layout.Column = 1;
            app.LastButton.Text = 'Last';

            % Create NextButton
            app.NextButton = uibutton(app.GridLayout3, 'push');
            app.NextButton.ButtonPushedFcn = createCallbackFcn(app, @NextButtonPushed, true);
            app.NextButton.Enable = 'off';
            app.NextButton.Layout.Row = 1;
            app.NextButton.Layout.Column = 3;
            app.NextButton.Text = 'Next';

            % Create EventsSlider
            app.EventsSlider = uislider(app.GridLayout3);
            app.EventsSlider.Limits = [1 100];
            app.EventsSlider.ValueChangedFcn = createCallbackFcn(app, @EventsSliderValueChanged, true);
            app.EventsSlider.Enable = 'off';
            app.EventsSlider.Layout.Row = 1;
            app.EventsSlider.Layout.Column = 2;
            app.EventsSlider.Value = 1;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout2);
            app.GridLayout4.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};
            app.GridLayout4.RowHeight = {'1x'};
            app.GridLayout4.Padding = [0 0 0 0];
            app.GridLayout4.Layout.Row = 1;
            app.GridLayout4.Layout.Column = 1;

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout4);
            app.GridLayout7.ColumnWidth = {'1x'};
            app.GridLayout7.Padding = [0 0 0 0];
            app.GridLayout7.Layout.Row = 1;
            app.GridLayout7.Layout.Column = 5;

            % Create SWThresholdLabel
            app.SWThresholdLabel = uilabel(app.GridLayout7);
            app.SWThresholdLabel.HorizontalAlignment = 'center';
            app.SWThresholdLabel.FontWeight = 'bold';
            app.SWThresholdLabel.Layout.Row = 2;
            app.SWThresholdLabel.Layout.Column = 1;
            app.SWThresholdLabel.Text = 'SW Threshold';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.GridLayout4);
            app.GridLayout8.ColumnWidth = {'1x'};
            app.GridLayout8.Padding = [0 0 0 0];
            app.GridLayout8.Layout.Row = 1;
            app.GridLayout8.Layout.Column = 4;

            % Create RSpinner
            app.RSpinner = uispinner(app.GridLayout8);
            app.RSpinner.Enable = 'off';
            app.RSpinner.Layout.Row = 1;
            app.RSpinner.Layout.Column = 1;

            % Create RThresholdLabel
            app.RThresholdLabel = uilabel(app.GridLayout8);
            app.RThresholdLabel.HorizontalAlignment = 'center';
            app.RThresholdLabel.FontWeight = 'bold';
            app.RThresholdLabel.Layout.Row = 2;
            app.RThresholdLabel.Layout.Column = 1;
            app.RThresholdLabel.Text = 'R Threshold';

            % Create SharpwaveCheckBox
            app.SharpwaveCheckBox = uicheckbox(app.GridLayout4);
            app.SharpwaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @SharpwaveCheckBoxValueChanged, true);
            app.SharpwaveCheckBox.Enable = 'off';
            app.SharpwaveCheckBox.Text = 'TP SW';
            app.SharpwaveCheckBox.FontWeight = 'bold';
            app.SharpwaveCheckBox.Layout.Row = 1;
            app.SharpwaveCheckBox.Layout.Column = 1;

            % Create RippleCheckBox
            app.RippleCheckBox = uicheckbox(app.GridLayout4);
            app.RippleCheckBox.ValueChangedFcn = createCallbackFcn(app, @RippleCheckBoxValueChanged, true);
            app.RippleCheckBox.Enable = 'off';
            app.RippleCheckBox.Text = 'TP R';
            app.RippleCheckBox.FontWeight = 'bold';
            app.RippleCheckBox.Layout.Row = 1;
            app.RippleCheckBox.Layout.Column = 2;

            % Create SaveButton
            app.SaveButton = uibutton(app.GridLayout4, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.WordWrap = 'on';
            app.SaveButton.Enable = 'off';
            app.SaveButton.Layout.Row = 1;
            app.SaveButton.Layout.Column = 3;
            app.SaveButton.Text = 'Save thresholds review';

            % Create SaveButton
            app.SaveButton = uibutton(app.GridLayout4, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.WordWrap = 'on';
            app.SaveButton.Enable = 'off';
            app.SaveButton.Layout.Row = 1;
            app.SaveButton.Layout.Column = 3;
            app.SaveButton.Text = 'Save thresholds review';

            % Create SaveButton
            app.SaveButton = uibutton(app.GridLayout4, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.WordWrap = 'on';
            app.SaveButton.Enable = 'off';
            app.SaveButton.Layout.Row = 1;
            app.SaveButton.Layout.Column = 3;
            app.SaveButton.Text = 'Save thresholds review';

            % Create RippleCheckBox
            app.RippleCheckBox = uicheckbox(app.GridLayout4);
            app.RippleCheckBox.ValueChangedFcn = createCallbackFcn(app, @RippleCheckBoxValueChanged, true);
            app.RippleCheckBox.Enable = 'off';
            app.RippleCheckBox.Text = 'TP R';
            app.RippleCheckBox.FontWeight = 'bold';
            app.RippleCheckBox.Layout.Row = 1;
            app.RippleCheckBox.Layout.Column = 2;

            % Create RippleCheckBox
            app.RippleCheckBox = uicheckbox(app.GridLayout4);
            app.RippleCheckBox.ValueChangedFcn = createCallbackFcn(app, @RippleCheckBoxValueChanged, true);
            app.RippleCheckBox.Enable = 'off';
            app.RippleCheckBox.Text = 'TP R';
            app.RippleCheckBox.FontWeight = 'bold';
            app.RippleCheckBox.Layout.Row = 1;
            app.RippleCheckBox.Layout.Column = 2;

            % Create SharpwaveCheckBox
            app.SharpwaveCheckBox = uicheckbox(app.GridLayout4);
            app.SharpwaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @SharpwaveCheckBoxValueChanged, true);
            app.SharpwaveCheckBox.Enable = 'off';
            app.SharpwaveCheckBox.Text = 'TP SW';
            app.SharpwaveCheckBox.FontWeight = 'bold';
            app.SharpwaveCheckBox.Layout.Row = 1;
            app.SharpwaveCheckBox.Layout.Column = 1;

            % Create SharpwaveCheckBox
            app.SharpwaveCheckBox = uicheckbox(app.GridLayout4);
            app.SharpwaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @SharpwaveCheckBoxValueChanged, true);
            app.SharpwaveCheckBox.Enable = 'off';
            app.SharpwaveCheckBox.Text = 'TP SW';
            app.SharpwaveCheckBox.FontWeight = 'bold';
            app.SharpwaveCheckBox.Layout.Row = 1;
            app.SharpwaveCheckBox.Layout.Column = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ThresholdReviewer_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end